"""Fused extract-and-infer engine: caches every gene tree's geometry once, then
reads each triplet's observations from it and runs inference in bounded
batches, serially or across workers.
"""

from io import StringIO
from typing import NamedTuple

from Bio import Phylo

from .config import DEFAULT_TREE_HEIGHT_CALCULATION_STRATEGY
from .inference import (
    _apply_triplet_result_p_value_correction,
    analyze_triplet_from_observations,
)
from .parallel import map_ordered, resolve_worker_count
from .triplet_geometry import (
    build_taxon_index,
    build_triplet_geometry,
    geometry_observation,
)

# How many triplets a worker holds observations for at once. One observation
# is a tuple and a float, about 96 bytes, so one triplet's list is ~260 KB at
# 2,700 gene trees and a whole 20,000-triplet chunk would be over 5 GB per
# worker. Extracting a batch at a time keeps that buffer in the tens of MB
# while the results, which are small, still accumulate for the chunk.
_EXTRACTION_BATCH_SIZE = 128

def _split_inference_kwargs(inference_kwargs):
    """Split the extraction-time keys from the per-triplet analysis kwargs.

    Args:
        inference_kwargs: Combined kwargs, optionally containing
            ``tree_height_calculation_strategy`` and
            ``collect_summary_statistics``.

    Returns:
        A tuple ``(strategy, collect_summary_statistics, analysis_kwargs)`` where
        ``analysis_kwargs`` is a copy without the extraction-time keys (safe to
        forward to :func:`analyze_triplet_from_observations`).
    """
    kwargs = dict(inference_kwargs or {})
    strategy = kwargs.pop(
        "tree_height_calculation_strategy", DEFAULT_TREE_HEIGHT_CALCULATION_STRATEGY
    )
    collect_summary_statistics = kwargs.pop("collect_summary_statistics", False)
    return strategy, collect_summary_statistics, kwargs


def _chunk_list(items, chunk_size):
    """Yield successive fixed-size chunks from a list.

    Args:
        items: The list to chunk.
        chunk_size: Maximum size of each chunk.

    Yields:
        Successive sublists of at most ``chunk_size`` items.
    """
    for i in range(0, len(items), chunk_size):
        yield items[i : i + chunk_size]


class RunGeometry(NamedTuple):
    """Every gene tree's cached geometry, built once for the whole run.

    Attributes:
        taxon_index: Taxon-to-position map spanning every triplet taxon, shared
            by all the cached geometries.
        geometries: One :class:`~.triplet_geometry.TripletGeometry` per gene
            tree, positionally aligned with the gene-tree list.
    """

    taxon_index: dict
    geometries: list


def _build_tree_geometry(taxon_index, newick_str):
    """Parse one processed gene tree and cache its geometry.

    Args:
        taxon_index: The run's shared taxon-to-position map.
        newick_str: The gene tree's processed Newick line.

    Returns:
        The tree's :class:`~.triplet_geometry.TripletGeometry`.
    """
    return build_triplet_geometry(
        Phylo.read(StringIO(newick_str), "newick"), taxon_index
    )


def build_run_geometry(gene_trees, triplets, processes=1):
    """Parse and cache every gene tree once, for reuse across all triplets.

    Each tree is read back from its processed Newick line, so the run measures
    exactly the branch lengths the processed-trees file records.

    Args:
        gene_trees: Gene-tree Newick strings.
        triplets: Every triplet the run will analyze, used to size the shared
            taxon index.
        processes: Worker processes; ``0`` means every CPU available to this
            process.

    Returns:
        A tuple ``(run_geometry, worker_cpu_seconds)``: the
        :class:`RunGeometry` holding the shared taxon index and one cached
        geometry per gene tree, and the CPU time spent in worker processes.
    """
    taxon_index = build_taxon_index(triplets)
    geometries, worker_cpu_seconds = map_ordered(
        _build_tree_geometry,
        gene_trees,
        resolve_worker_count(len(gene_trees), processes),
        shared=taxon_index,
    )
    return RunGeometry(taxon_index, geometries), worker_cpu_seconds


def _extract_chunk_observations(
    triplet_chunk,
    run_geometry,
    strategy,
    collect_summary_statistics=False,
):
    """Read every triplet's observations out of the cached gene-tree geometries.

    Args:
        triplet_chunk: Iterable of triplets in this chunk.
        run_geometry: The :class:`RunGeometry` cached for the whole run.
        strategy: Tree-height strategy used to compute each observation.
        collect_summary_statistics: When ``True``, each observation also carries
            its per-tree summary metrics.

    Returns:
        A dict mapping each triplet to its list of
        ``(topology, tree_height, summary_metrics)`` observations, in gene-tree
        order.
    """
    observations = {triplet: [] for triplet in triplet_chunk}
    taxon_index = run_geometry.taxon_index
    positions = {
        triplet: tuple(taxon_index[label] for label in triplet)
        for triplet in triplet_chunk
    }

    for geometry in run_geometry.geometries:
        for triplet in triplet_chunk:
            observation = geometry_observation(
                geometry,
                positions[triplet],
                strategy,
                collect_summary_statistics,
            )
            if observation is not None:
                observations[triplet].append(observation)
    return observations


def _analyze_chunk(
    triplet_chunk,
    run_geometry,
    species_triplet_trees,
    strategy,
    collect_summary_statistics,
    analysis_kwargs,
):
    """Run the fused extract-then-infer loop for one triplet chunk.

    Args:
        triplet_chunk: Iterable of triplets in this chunk.
        run_geometry: The run-wide geometry cache.
        species_triplet_trees: Mapping of triplet to its species subtree Newick.
        strategy: Tree-height strategy used to compute observations.
        collect_summary_statistics: When ``True``, gather per-triplet summary
            statistics during extraction.
        analysis_kwargs: Keyword arguments forwarded to
            :func:`analyze_triplet_from_observations`.

    Returns:
        A list of ``TripletPipelineResult`` objects for the chunk.
    """
    results = []
    for batch in _chunk_list(list(triplet_chunk), _EXTRACTION_BATCH_SIZE):
        observations = _extract_chunk_observations(
            batch,
            run_geometry,
            strategy,
            collect_summary_statistics,
        )
        for triplet in batch:
            results.append(
                analyze_triplet_from_observations(
                    triplet,
                    observations[triplet],
                    species_subtree=species_triplet_trees.get(triplet),
                    **analysis_kwargs,
                )
            )
    return results


def _analyze_chunk_shared(shared, triplet_chunk):
    """Run the fused loop for one triplet chunk against the shared run state.

    Args:
        shared: ``(run_geometry, species_triplet_trees, inference_kwargs)``.
        triplet_chunk: Iterable of triplets in this chunk.

    Returns:
        A list of ``TripletPipelineResult`` objects for the chunk.
    """
    run_geometry, species_triplet_trees, inference_kwargs = shared
    strategy, collect_summary, analysis_kwargs = _split_inference_kwargs(
        inference_kwargs
    )
    return _analyze_chunk(
        triplet_chunk,
        run_geometry,
        species_triplet_trees,
        strategy,
        collect_summary,
        analysis_kwargs,
    )


def _run_taxon_mode(
    triplets, species_triplet_trees, run_geometry, worker_count, inference_kwargs
):
    """Dispatch triplet chunks across workers, each running the fused loop.

    Args:
        triplets: List of triplets to analyze.
        species_triplet_trees: Mapping of triplet to its species subtree Newick.
        run_geometry: The run-wide geometry cache, shared read-only to workers.
        worker_count: Number of worker processes.
        inference_kwargs: Combined inference kwargs forwarded to workers.

    Returns:
        A tuple ``(results, worker_cpu_seconds)`` where ``results`` is a flat
        list of ``TripletPipelineResult`` objects and ``worker_cpu_seconds`` is
        the total CPU time spent across all workers.
    """
    chunksize = max(1, len(triplets) // (worker_count * 4))
    chunk_results, worker_cpu_seconds = map_ordered(
        _analyze_chunk_shared,
        _chunk_list(triplets, chunksize),
        worker_count,
        shared=(run_geometry, species_triplet_trees, inference_kwargs),
        chunksize=1,
    )
    results = [result for chunk in chunk_results for result in chunk]
    return results, worker_cpu_seconds


def stream_triplet_results(
    triplets,
    species_triplet_trees,
    gene_trees,
    *,
    processes=0,
    inference_kwargs=None,
    p_value_correction="no",
    return_worker_cpu=False,
):
    """Fuse extraction and inference over all triplets and apply run-wide correction.

    Resolves the worker count, caches every gene tree's geometry once, runs the
    fused engine over triplet chunks, then applies the run-wide p-value
    correction in a single pass.

    Args:
        triplets: List of triplets to analyze.
        species_triplet_trees: Mapping of triplet to its species subtree Newick.
        gene_trees: Gene-tree Newick strings.
        processes: Requested worker count; ``0`` means every CPU available to
            this process.
        inference_kwargs: Keyword arguments for the per-triplet analysis (also
            the source of ``tree_height_calculation_strategy`` used at
            extraction and ``alpha_dct``/``alpha_ks`` used for correction).
        p_value_correction: Correction method applied across all triplets.
        return_worker_cpu: When ``True``, also return the CPU seconds spent in
            pool workers (which the caller's own ``process_time`` cannot see).

    Returns:
        A list of corrected ``TripletPipelineResult`` objects (empty if there
        are no triplets). When ``return_worker_cpu`` is ``True``, a tuple
        ``(results, worker_cpu_seconds)`` instead.
    """
    inference_kwargs = dict(inference_kwargs or {})
    if not triplets:
        return ([], 0.0) if return_worker_cpu else []

    # Bootstrap iterations need the same correction the point estimate gets, so
    # the method and the family size travel with the per-triplet analysis.
    inference_kwargs["p_value_correction"] = p_value_correction
    inference_kwargs["family_size"] = len(triplets)

    worker_count = resolve_worker_count(len(triplets), processes)

    # Cache every gene tree's geometry once for the whole run. Chunks reuse it,
    # so a tree is parsed once rather than once per chunk per worker.
    strategy, collect_summary, analysis_kwargs = _split_inference_kwargs(
        inference_kwargs
    )
    run_geometry, worker_cpu_seconds = build_run_geometry(
        gene_trees, triplets, processes=processes
    )

    if worker_count <= 1:
        # Serial: work runs in the parent, so its CPU is captured by the caller.
        results = _analyze_chunk(
            triplets,
            run_geometry,
            species_triplet_trees,
            strategy,
            collect_summary,
            analysis_kwargs,
        )
    else:
        results, inference_worker_cpu = _run_taxon_mode(
            triplets,
            species_triplet_trees,
            run_geometry,
            worker_count,
            inference_kwargs,
        )
        worker_cpu_seconds += inference_worker_cpu

    corrected = _apply_triplet_result_p_value_correction(
        results,
        alpha_dct=inference_kwargs.get("alpha_dct"),
        alpha_ks=inference_kwargs.get("alpha_ks"),
        method=p_value_correction,
    )

    if return_worker_cpu:
        return corrected, worker_cpu_seconds
    return corrected
