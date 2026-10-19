"""Simulate a species tree and gene trees with msprime for benchmarking GhostParser.

    python sim_script/simulate_trees.py -c sim_script/configs/basic_newick.yaml

See sim_script/README.md for every config key.
"""

from __future__ import annotations

import argparse
import csv
import shutil
import sys
import traceback
from pathlib import Path

import numpy as np
import yaml

from sim_config import SimConfig, SimError, SimInputError, load_config

EXIT_OK, EXIT_INTERNAL, EXIT_INTERRUPTED = 0, 70, 130


def _require_packages() -> None:
    missing = []
    for module in ("msprime", "tskit", "Bio", "matplotlib"):
        try:
            __import__(module)
        except ImportError:
            missing.append(module)
    if missing:
        raise SimInputError(
            f"Missing package(s): {', '.join(missing)}; install them with "
            "pip install -r sim_script/requirements.txt"
        )


def prepare_output_dir(path: Path, overwrite: bool) -> Path:
    """Reset ``path`` when ``overwrite`` is set; otherwise use the first free ``<name>_<n>`` sibling."""
    if path.exists() and overwrite:
        shutil.rmtree(path) if path.is_dir() else path.unlink()
    if path.exists():
        n = 1
        while path.with_name(f"{path.name}_{n}").exists():
            n += 1
        path = path.with_name(f"{path.name}_{n}")
    path.mkdir(parents=True)
    return path


def _ancestry_model(spec):
    import msprime

    if isinstance(spec, str):
        return spec
    (name, params), = spec.items()
    if name == "beta":
        return msprime.BetaCoalescent(**params)
    return msprime.DiracCoalescent(**params)


def _length_factor(cfg: SimConfig, origin_age: float) -> float:
    """Output length per generation, at the mean clock rate for the substitution units.

    Args:
        cfg: The validated config.
        origin_age: Age of the species tree's root in generations; ``normalized``
            divides by it, so the species tree's origin-to-present depth is 1.

    Returns:
        The factor, ``output.scale`` included.
    """
    units = cfg["output"]["branch_length_units"]
    scale = cfg["output"]["scale"]
    if units == "normalized":
        return scale / origin_age
    if units == "substitutions":
        return scale * cfg["clock"]["mean_rate"]
    if units == "time_units":
        return scale / cfg.generations_per_unit
    if units == "generations":
        return scale
    pop = cfg["population"]
    return scale / (pop["ploidy"] * pop["size"])


def _gene_tree_newick(tree, labels: dict[int, str], factor: float, clock: dict | None,
                      rng: np.random.Generator, decimals: int) -> str:
    """Write one tskit tree as Newick, applying the clock when ``clock`` is given.

    Each branch's rate is its own lognormal draw (an uncorrelated relaxed clock)
    with mean 1 relative to the mean rate, so the factor carries the units.
    """
    locus = 1.0
    if clock is not None and clock["among_loci_gamma_shape"] is not None:
        shape = clock["among_loci_gamma_shape"]
        locus = rng.gamma(shape, 1.0 / shape)
    sigma = clock["sigma"] if clock is not None and clock["model"] == "ucln" else 0.0

    def write(u: int) -> str:
        kids = tree.children(u)
        body = labels[u] if not kids else "(" + ",".join(write(c) for c in kids) + ")"
        if u == tree.root:
            return body
        length = tree.branch_length(u) * factor * locus
        if sigma > 0:
            length *= rng.lognormal(-sigma * sigma / 2, sigma)
        return f"{body}:{length:.{decimals}f}"

    return write(tree.root) + ";"


def _write_gene_trees(cfg: SimConfig, demography, root, path: Path, ancestry_seed: int,
                      rate_rng: np.random.Generator) -> None:
    import msprime

    from demography import population_name
    from sim_config import CLOCK_UNITS
    from species_tree import find_tip

    sampled = list(cfg.species.sampled)
    pop_of = {population_name(find_tip(root, name)): name for name in sampled}
    samples = [msprime.SampleSet(1, population=pop, ploidy=1) for pop in pop_of]
    locus = cfg["locus"]
    replicates = msprime.sim_ancestry(
        samples=samples, demography=demography, sequence_length=locus["sequence_length"],
        recombination_rate=locus["recombination_rate"], ploidy=cfg["population"]["ploidy"],
        model=_ancestry_model(cfg["ancestry_model"]), num_replicates=cfg["gene_tree_count"],
        random_seed=ancestry_seed,
    )
    factor = _length_factor(cfg, root.age)
    clock = cfg["clock"] if cfg["output"]["branch_length_units"] in CLOCK_UNITS else None
    decimals = cfg["output"]["decimals"]
    with path.open("w") as handle:
        for ts in replicates:
            # One sample per species: label each by the population it was drawn from.
            labels = {u: pop_of[ts.population(ts.node(u).population).metadata["name"]] for u in ts.samples()}
            if locus["tree_choice"] == "first":
                tree = ts.first()
            else:
                tree = ts.at(ts.sequence_length / 2)
            handle.write(_gene_tree_newick(tree, labels, factor, clock, rate_rng, decimals) + "\n")


def _write_species_trees(cfg: SimConfig, root, rep_dir: Path) -> None:
    from species_tree import prune, to_newick

    factor = _length_factor(cfg, root.age)
    decimals = cfg["output"]["decimals"]

    def length(node) -> float:
        return (node.parent.age - node.age) * factor

    (rep_dir / "species_full.tree").write_text(to_newick(root, length, decimals) + "\n")
    sampled = prune(root, set(cfg.species.sampled))
    (rep_dir / "species.tree").write_text(to_newick(sampled, length, decimals) + "\n")


def _write_events(cfg: SimConfig, events, path: Path) -> None:
    from demography import EVENT_COLUMNS

    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=EVENT_COLUMNS, delimiter="\t")
        writer.writeheader()
        for event in events:
            writer.writerow(event.row(cfg.generations_per_unit))


def _seed_int(seq: np.random.SeedSequence) -> int:
    """A seed msprime accepts (1 to 2**32 - 1) drawn from a seed sequence."""
    return int(seq.generate_state(1, dtype=np.uint32)[0]) % (2**32 - 1) + 1


def run(cfg: SimConfig) -> Path:
    """Simulate every replicate into the output directory.

    Args:
        cfg: The validated config.

    Returns:
        The output directory used.

    Raises:
        SimError: A replicate's tree or events cannot be built from the config.
    """
    from demography import build_demography, resolve_events
    from species_tree import build_species_tree

    master = np.random.SeedSequence(cfg["seed"])
    seed = master.entropy
    # Build every replicate's scenario before touching the disk, so a config
    # that fails on some replicate leaves no half-written output behind.
    scenarios = []
    for rep_seq in master.spawn(cfg["replicates"]):
        tree_seq, event_seq, ancestry_seq, rate_seq = rep_seq.spawn(4)
        root = build_species_tree(cfg, np.random.default_rng(tree_seq))
        events = resolve_events(cfg, root, np.random.default_rng(event_seq))
        scenarios.append((root, events, build_demography(cfg, root, events), ancestry_seq, rate_seq))

    out = prepare_output_dir(cfg.output_dir, cfg["overwrite"])
    for i, (root, events, demography, ancestry_seq, rate_seq) in enumerate(scenarios):
        rep_dir = out / f"rep_{i:03d}"
        rep_dir.mkdir()
        _write_species_trees(cfg, root, rep_dir)
        _write_events(cfg, events, rep_dir / "events.tsv")
        with (rep_dir / "demography.txt").open("w") as handle:
            demography.debug().print_history(output=handle)
        _write_gene_trees(cfg, demography, root, rep_dir / "genes.tree", _seed_int(ancestry_seq),
                          np.random.default_rng(rate_seq))

        record = dict(cfg.raw, seed=seed, output_dir=str(out))
        (rep_dir / "resolved_config.yaml").write_text(yaml.safe_dump(record, sort_keys=False))
        if cfg["output"]["write_ghostparser_config"]:
            gp = {
                "species_tree_path": str(rep_dir / "species.tree"),
                "gene_trees_path": str(rep_dir / "genes.tree"),
                "outgroup": list(cfg.species.outgroup),
                "output_folder": str(rep_dir / "ghostparser_results"),
            }
            (rep_dir / "ghostparser_config.yaml").write_text(yaml.safe_dump(gp, sort_keys=False))
        if cfg["output"]["plot"]:
            _plot(cfg, root, events, rep_dir, i, seed)
        print(f"rep_{i:03d}: {len(events)} event(s), {cfg['gene_tree_count']} gene trees -> {rep_dir}")
    return out


def _plot(cfg: SimConfig, root, events, rep_dir: Path, replicate: int, seed: int) -> None:
    from plot_events import plot_scenario

    tree = cfg.tree
    model = tree["model"]
    if model == "birth_death":
        model += f" (birth {tree['birth_rate']:g}, death {tree['death_rate']:g})"
    locus = cfg["locus"]
    caption = (f"species tree: {model} | seed {seed} | Ne {cfg['population']['size']:g} | "
               f"{cfg['gene_tree_count']} gene trees | locus {locus['sequence_length']:g} bp, "
               f"recombination {locus['recombination_rate']:g} | {cfg['ancestry_model']} | "
               f"clock {cfg['clock']['model']}")
    path = rep_dir / f"scenario.{cfg['output']['plot']}"
    try:
        plot_scenario(cfg, root, events, path, f"Replicate {replicate}", caption)
    except Exception as exc:  # a figure must never cost the simulated data
        print(f"Warning: could not draw {path.name}: {exc}", file=sys.stderr)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("-c", "--config", required=True, help="YAML or JSON config file")
    parser.add_argument("-o", "--output-dir", help="Output directory (overrides output_dir)")
    parser.add_argument("--seed", type=int, help="Master seed (overrides seed)")
    parser.add_argument("--overwrite", action="store_true", default=None,
                        help="Reset an existing output directory instead of writing beside it")
    parser.add_argument("--debug", action="store_true", help="Show the full traceback on failure")
    args = parser.parse_args(argv)
    try:
        _require_packages()
        cfg = load_config(args.config, {"output_dir": args.output_dir, "seed": args.seed,
                                        "overwrite": args.overwrite})
        out = run(cfg)
    except KeyboardInterrupt:
        if args.debug:
            raise
        print("Interrupted.", file=sys.stderr)
        return EXIT_INTERRUPTED
    except SimError as exc:
        if args.debug:
            raise
        label = "Config error" if exc.exit_status == 2 else "Error"
        print(f"{label}: {exc}", file=sys.stderr)
        return exc.exit_status
    except Exception:
        if args.debug:
            raise
        print("Internal error (rerun with --debug for the traceback):", file=sys.stderr)
        print(traceback.format_exc().strip().splitlines()[-1], file=sys.stderr)
        return EXIT_INTERNAL
    print(f"Done: {out}")
    return EXIT_OK


if __name__ == "__main__":
    sys.exit(main())
