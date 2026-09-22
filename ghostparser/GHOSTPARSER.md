# GhostParser Package Guide

The `ghostparser` package is the introgression engine plus an optional
machine-learning subpackage, sharing a thin configuration layer.

```
ghostparser/
  orchestrator/     the introgression engine and its consolidation stage
  ml/               optional multi-label classifiers over summary_statistics.tsv
  config.py         configuration helpers shared by every module
  cli_config.py     config file + command line resolution shared by every module
  triplet_utils.py  triplet topology helpers
  __main__.py       usage banner for `python -m ghostparser`
```

| Module | Purpose | Guide |
| --- | --- | --- |
| `ghostparser.orchestrator` | Roots the trees, reads every species triplet out of every gene tree, runs the three-gate test cascade, and draws the introgression maps. The primary entry point. | [orchestrator/ORCHESTRATOR.md](orchestrator/ORCHESTRATOR.md) |
| `ghostparser.ml` | Trains a random forest or a multi-label KNN on a run's `summary_statistics.tsv` and tunes their hyperparameters. Needs `pip install .[ml]`; Weights & Biases logging is the `wandb` extra. | [ml/ML.md](ml/ML.md) |

The runnable commands are `python -m ghostparser.orchestrator`,
`python -m ghostparser.ml.random_forest`, `python -m ghostparser.ml.multi_knn`
and `python -m ghostparser.ml.hyper_tune`; `python -m ghostparser` and
`python -m ghostparser.ml` print usage.

The two modules own their own defaults, choices and validation, so their
settings can diverge; what they share is the configuration error type, path
resolution, config-file loading, the overwrite flag and output-directory
preparation, and one resolution rule: a config file supplies the settings and
any flag given beside it overrides the file's value. Every key is documented
in [CONFIG.md](../CONFIG.md).
