# Changelog

## v0.1.0 - Oct 8, 2026

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
  simulated training statistics, stored with Git LFS.
