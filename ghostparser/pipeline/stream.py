"""Fused extract-and-infer streaming engine for the pipeline.

Extracts a chunk of triplets' observations in one gene-tree parse pass, runs
inference immediately, and never materializes the global intermediate
structure. Observations (topology + tree height) are computed directly from the
extracted subtree objects, so no serialize-then-reparse round trip happens. Also
resolves the parallelization mode and applies the run-wide p-value correction.
"""

from __future__ import annotations

import time
from multiprocessing import cpu_count

import dendropy

from ghostparser.config import DEFAULT_TREE_HEIGHT_CALCULATION_STRATEGY

from .config import AUTO_GENE_TREES_THRESHOLD, AUTO_TAXA_SMALL_THRESHOLD
from .inference import (
    _apply_triplet_result_p_value_correction,
    analyze_triplet_from_observations,
    observation_from_subtree,
)
from .trees import _get_mp_context, extract_triplet_subtree

# Worker-global state shared read-only via fork/forkserver initializer.
_GENE_TREES: list[str] = []
_SPECIES_TRIPLET_TREES: dict[tuple[str, str, str], str] = {}
_STRATEGY: str = DEFAULT_TREE_HEIGHT_CALCULATION_STRATEGY
_ANALYSIS_KWARGS: dict = {}


def _split_inference_kwargs(inference_kwargs):
    """Split the tree-height strategy (used at extraction) from analysis kwargs.

    Args:
        inference_kwargs: Combined kwargs, optionally containing
            ``tree_height_calculation_strategy``.

    Returns:
        A tuple ``(strategy, analysis_kwargs)`` where ``analysis_kwargs`` is a
        copy without the strategy key (safe to forward to
        :func:`analyze_triplet_from_observations`).
    """
    kwargs = dict(inference_kwargs or {})
    strategy = kwargs.pop(
        "tree_height_calculation_strategy", DEFAULT_TREE_HEIGHT_CALCULATION_STRATEGY
    )
    return strategy, kwargs


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


def _calculate_worker_count(total_items, processes):
    """Resolve the worker count from the requested processes and item count.

    Args:
        total_items: Number of work items available.
        processes: Requested worker count; ``0`` means all cores.

    Returns:
        The clamped worker count (at least 1, at most ``total_items``).
    """
    available_cpus = cpu_count()
    worker_count = processes or available_cpus
    return max(1, min(worker_count, total_items))


def resolve_parallelization_mode(mode, n_taxa, n_gene_trees):
    """Resolve the ``auto`` mode to ``gene`` or ``taxon`` from input size.

    Args:
        mode: The requested mode (``auto``, ``taxon``, or ``gene``).
        n_taxa: Number of distinct ingroup taxa.
        n_gene_trees: Number of gene trees.

    Returns:
        ``taxon`` or ``gene`` (``mode`` unchanged when it is not ``auto``).
    """
    if mode != "auto":
        return mode
    if n_taxa < AUTO_TAXA_SMALL_THRESHOLD or n_gene_trees > AUTO_GENE_TREES_THRESHOLD:
        return "gene"
    return "taxon"


def _extract_chunk_observations(triplet_chunk, gene_trees, strategy):
    """Compute every triplet's observations in one pass over the gene trees.

    Args:
        triplet_chunk: Iterable of triplets in this chunk.
        gene_trees: Iterable of gene-tree Newick strings.
        strategy: Tree-height strategy used to compute each observation.

    Returns:
        A dict mapping each triplet to its list of ``(topology, tree_height)``
        observations.
    """
    observations = {triplet: [] for triplet in triplet_chunk}
    for newick_str in gene_trees:
        tree = dendropy.Tree.get(
            data=newick_str, schema="newick", preserve_underscores=True
        )
        tree_taxa = {taxon.label for taxon in tree.taxon_namespace if taxon.label}
        for triplet in triplet_chunk:
            if not set(triplet).issubset(tree_taxa):
                continue
            subtree = extract_triplet_subtree(tree, triplet)
            if subtree is None:
                continue
            observation = observation_from_subtree(subtree, triplet, strategy)
            if observation is not None:
                observations[triplet].append(observation)
    return observations


def _analyze_chunk(
    triplet_chunk, gene_trees, species_triplet_trees, strategy, analysis_kwargs
):
    """Run the fused extract-then-infer loop for one triplet chunk.

    Args:
        triplet_chunk: Iterable of triplets in this chunk.
        gene_trees: Iterable of gene-tree Newick strings.
        species_triplet_trees: Mapping of triplet to its species subtree Newick.
        strategy: Tree-height strategy used to compute observations.
        analysis_kwargs: Keyword arguments forwarded to
            :func:`analyze_triplet_from_observations`.

    Returns:
        A list of ``TripletPipelineResult`` objects for the chunk.
    """
    observations = _extract_chunk_observations(triplet_chunk, gene_trees, strategy)
    results = []
    for triplet in triplet_chunk:
        results.append(
            analyze_triplet_from_observations(
                triplet,
                observations[triplet],
                species_subtree=species_triplet_trees.get(triplet),
                **analysis_kwargs,
            )
        )
    return results


def _init_stream_worker(gene_trees, species_triplet_trees, inference_kwargs):
    """Seed worker globals with shared read-only state.

    Args:
        gene_trees: Gene-tree Newick strings, or ``None`` to keep the value the
            parent already set via fork.
        species_triplet_trees: Mapping of triplet to its species subtree Newick.
        inference_kwargs: Combined inference kwargs (split into strategy and
            analysis kwargs).
    """
    global _GENE_TREES, _SPECIES_TRIPLET_TREES, _STRATEGY, _ANALYSIS_KWARGS
    if gene_trees is not None:
        _GENE_TREES = gene_trees
    _SPECIES_TRIPLET_TREES = species_triplet_trees or {}
    _STRATEGY, _ANALYSIS_KWARGS = _split_inference_kwargs(inference_kwargs)


def _analyze_chunk_worker(triplet_chunk):
    """Analyze one triplet chunk against the worker-global gene trees.

    Args:
        triplet_chunk: Iterable of triplets in this chunk.

    Returns:
        A tuple ``(results, worker_cpu_seconds)`` where ``results`` is a list of
        ``TripletPipelineResult`` objects and ``worker_cpu_seconds`` is the CPU
        time this worker spent on the chunk.
    """
    start_cpu = time.process_time()
    results = _analyze_chunk(
        triplet_chunk, _GENE_TREES, _SPECIES_TRIPLET_TREES, _STRATEGY, _ANALYSIS_KWARGS
    )
    return results, time.process_time() - start_cpu


def _extract_triplet_range_worker(args):
    """Compute one triplet's observations from a gene-tree index range.

    Args:
        args: A tuple ``(triplet, lo, hi)`` selecting ``_GENE_TREES[lo:hi]``.

    Returns:
        A tuple ``(observations, worker_cpu_seconds)`` where ``observations`` is
        a list of ``(topology, tree_height)`` tuples and ``worker_cpu_seconds``
        is the CPU time this worker spent on the range.
    """
    start_cpu = time.process_time()
    triplet, lo, hi = args
    out = []
    triplet_set = set(triplet)
    for newick_str in _GENE_TREES[lo:hi]:
        tree = dendropy.Tree.get(
            data=newick_str, schema="newick", preserve_underscores=True
        )
        tree_taxa = {taxon.label for taxon in tree.taxon_namespace if taxon.label}
        if not triplet_set.issubset(tree_taxa):
            continue
        subtree = extract_triplet_subtree(tree, triplet)
        if subtree is None:
            continue
        observation = observation_from_subtree(subtree, triplet, _STRATEGY)
        if observation is not None:
            out.append(observation)
    return out, time.process_time() - start_cpu


def _run_taxon_mode(triplets, species_triplet_trees, gene_trees, worker_count, inference_kwargs):
    """Dispatch triplet chunks across workers, each running the fused loop.

    Args:
        triplets: List of triplets to analyze.
        species_triplet_trees: Mapping of triplet to its species subtree Newick.
        gene_trees: Gene-tree Newick strings, shared read-only to workers.
        worker_count: Number of worker processes.
        inference_kwargs: Combined inference kwargs forwarded to workers.

    Returns:
        A tuple ``(results, worker_cpu_seconds)`` where ``results`` is a flat
        list of ``TripletPipelineResult`` objects and ``worker_cpu_seconds`` is
        the total CPU time spent across all workers.
    """
    chunksize = max(1, len(triplets) // (worker_count * 4))
    triplet_chunks = list(_chunk_list(triplets, chunksize))

    ctx = _get_mp_context()
    init_gene_trees = gene_trees
    if hasattr(ctx, "get_start_method") and ctx.get_start_method() == "fork":
        global _GENE_TREES
        _GENE_TREES = gene_trees
        init_gene_trees = None

    with ctx.Pool(
        processes=worker_count,
        initializer=_init_stream_worker,
        initargs=(init_gene_trees, species_triplet_trees, inference_kwargs),
    ) as pool:
        chunk_payloads = list(
            pool.imap(_analyze_chunk_worker, iter(triplet_chunks), chunksize=1)
        )

    results = [result for payload in chunk_payloads for result in payload[0]]
    worker_cpu_seconds = sum(payload[1] for payload in chunk_payloads)
    return results, worker_cpu_seconds


def _run_gene_mode(triplets, species_triplet_trees, gene_trees, worker_count, inference_kwargs):
    """Process triplets serially, parallelizing per-gene-tree extraction within each.

    Args:
        triplets: List of triplets to analyze.
        species_triplet_trees: Mapping of triplet to its species subtree Newick.
        gene_trees: Gene-tree Newick strings, shared read-only to workers.
        worker_count: Number of worker processes.
        inference_kwargs: Combined inference kwargs (strategy used by extraction
            workers, analysis kwargs used by the parent).

    Returns:
        A tuple ``(results, worker_cpu_seconds)`` where ``results`` is a list of
        ``TripletPipelineResult`` objects in triplet order and
        ``worker_cpu_seconds`` is the total CPU time spent across all extraction
        workers. Per-triplet inference and bootstrap run in the parent, so their
        CPU is captured by the caller's own timing.
    """
    _strategy, analysis_kwargs = _split_inference_kwargs(inference_kwargs)
    n_gene_trees = len(gene_trees)
    range_chunk = max(1, n_gene_trees // (worker_count * 4))
    ranges = [
        (lo, min(lo + range_chunk, n_gene_trees))
        for lo in range(0, n_gene_trees, range_chunk)
    ]

    ctx = _get_mp_context()
    init_gene_trees = gene_trees
    if hasattr(ctx, "get_start_method") and ctx.get_start_method() == "fork":
        global _GENE_TREES
        _GENE_TREES = gene_trees
        init_gene_trees = None

    results = []
    worker_cpu_seconds = 0.0
    with ctx.Pool(
        processes=worker_count,
        initializer=_init_stream_worker,
        initargs=(init_gene_trees, species_triplet_trees, inference_kwargs),
    ) as pool:
        for triplet in triplets:
            tasks = [(triplet, lo, hi) for lo, hi in ranges]
            payloads = pool.map(_extract_triplet_range_worker, tasks)
            observations = [obs for payload in payloads for obs in payload[0]]
            worker_cpu_seconds += sum(payload[1] for payload in payloads)
            results.append(
                analyze_triplet_from_observations(
                    triplet,
                    observations,
                    species_subtree=species_triplet_trees.get(triplet),
                    **analysis_kwargs,
                )
            )

    return results, worker_cpu_seconds


def stream_triplet_results(
    triplets,
    species_triplet_trees,
    gene_trees,
    *,
    mode="auto",
    processes=0,
    inference_kwargs=None,
    p_value_correction="no",
    return_worker_cpu=False,
):
    """Fuse extraction and inference over all triplets and apply run-wide correction.

    Resolves the mode and worker count, runs the fused engine (serial, ``taxon``,
    or ``gene``), then applies the run-wide p-value correction in a single pass.

    Args:
        triplets: List of triplets to analyze.
        species_triplet_trees: Mapping of triplet to its species subtree Newick.
        gene_trees: Gene-tree Newick strings.
        mode: Parallelization mode (``auto``, ``taxon``, or ``gene``).
        processes: Requested worker count; ``0`` means all cores.
        inference_kwargs: Keyword arguments for the per-triplet analysis (also
            the source of ``tree_height_calculation_strategy`` used at
            extraction and ``alpha_dct``/``alpha_ks``/``stats_backend`` used for
            correction).
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

    n_taxa = len({taxon for triplet in triplets for taxon in triplet})
    resolved_mode = resolve_parallelization_mode(mode, n_taxa, len(gene_trees))
    worker_count = _calculate_worker_count(len(triplets), processes)

    worker_cpu_seconds = 0.0
    if worker_count <= 1:
        # Serial: work runs in the parent, so its CPU is captured by the caller.
        strategy, analysis_kwargs = _split_inference_kwargs(inference_kwargs)
        results = _analyze_chunk(
            triplets, gene_trees, species_triplet_trees, strategy, analysis_kwargs
        )
    elif resolved_mode == "gene":
        results, worker_cpu_seconds = _run_gene_mode(
            triplets, species_triplet_trees, gene_trees, worker_count, inference_kwargs
        )
    else:
        results, worker_cpu_seconds = _run_taxon_mode(
            triplets, species_triplet_trees, gene_trees, worker_count, inference_kwargs
        )

    corrected = _apply_triplet_result_p_value_correction(
        results,
        alpha_dct=inference_kwargs.get("alpha_dct"),
        alpha_ks=inference_kwargs.get("alpha_ks"),
        method=p_value_correction,
        stats_backend=inference_kwargs.get("stats_backend", "standard"),
    )

    if return_worker_cpu:
        return corrected, worker_cpu_seconds
    return corrected
