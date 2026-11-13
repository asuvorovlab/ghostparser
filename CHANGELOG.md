# Changelog

## v0.1.2

- `outflow_introgression` is renamed `outflow_or_bidirectional_introgression`,
  and its `inference` reads `introgression from B to C with probable C to B`.
- Consolidation weights only the outflow edge and marks its reverse
  cell with a dashed inner box in both figures, black or white for contrast
  with the cell's fill; `introgression_matrix_sampled_probable_reverse.tsv`
  counts the triplets behind each mark.
- The sampled matrices are renamed from
  `introgression_matrix_inflow_outflow*.tsv` to
  `introgression_matrix_sampled*.tsv`.
- New `bootstrap_threshold` key (`-bt`, `--bootstrap-threshold`), off by
  default: the introgression maps take only triplets whose `bootstrap_value`
  clears their class's threshold, given as one value, four in the order no,
  inflow, outflow, ghost, or a mapping of `MSC`, `MSCi_in`, `MSCi_out` and
  `MSCig`. The results TSV keeps every triplet.

## v0.1.1 - Oct 30, 2026

- Consolidation writes two figures in place of `introgression_combined.png`.
  `introgression_bootstrap.png` keeps the bootstrap heatmap and draws every
  ghost target as an equal-length bar filled on the heatmap's colour scale, so
  the bar's colour rather than its length carries the support.
  `introgression_counts.png` shows supporting triplet counts: black heatmap
  cells labelled with their count and no colour scale, and ghost bars as long
  as their count. In both, a ghost bar's outline is dashed for ghost
  introgression only and solid when the taxon is also the target of a sampled
  edge.

## v0.1.0 - Oct 12, 2026

The initial release under this repository. Earlier versions were developed in
a private personal repository under the `asif256000` account, and the scripts
released with the 2025 bioRxiv preprint are archived in the
[legacy-2025 release](https://github.com/asuvorovlab/ghostparser/releases/tag/legacy-2025).

- A pip-installable Python package (Python 3.10 or later) that replaces the
  preprint's Python and R scripts, with no R or newick_utils dependency.
- `python -m ghostparser.orchestrator` roots the species tree and the gene
  trees on the outgroups, reads every ingroup triplet out of every gene tree,
  and classifies each triplet as no, inflow, outflow or ghost introgression,
  or ambiguous.
- Three tests per triplet: a discordant count test, a Kolmogorov-Smirnov
  tree-height test, and a studentized permutation test with a TOST
  equivalence step for the direction; the p-values are corrected across the
  run, and each call gets a bootstrap support value.
- The default settings match those used in the IDP preprint, and
  every one is configurable.
- Consolidated introgression maps (a heatmap of sampled introgression beside
  a bar chart of ghost targets, with their TSV matrices).
- Settings from command-line flags or a YAML or JSON config file, an input
  check that runs before any analysis, triplet and species filters, and
  one-line errors with distinct exit statuses. Runs in parallel on one
  machine, with identical results at any worker count under a fixed seed.
- The optional `ghostparser.ml` subpackage trains random-forest and
  multi-label KNN classifiers on a run's summary statistics and tunes their
  hyperparameters, optionally logging to Weights & Biases.
- The empirical datasets (Drosophila, Heliconius, Jaltomata, Thuja) and the
  simulated training statistics, kept in the repository rather than in the
  release downloads.
