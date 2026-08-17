import json

from ghostparser.ml.config import load_ml_config


def test_load_ml_config_defaults_target_column_to_class(tmp_path):
    config_path = tmp_path / "ml_config.json"
    config_path.write_text(
        json.dumps(
            {
                "input_path": "./results/summary_statistics.tsv",
                "output_dir": "./results/ml_out",
            }
        )
    )

    config = load_ml_config(str(config_path))

    assert config["target_column"] == "class"
    assert config["overwrite"] is True


def test_load_ml_config_accepts_explicit_class_target_column(tmp_path):
    config_path = tmp_path / "ml_config_explicit.json"
    config_path.write_text(
        json.dumps(
            {
                "input_path": "./results/summary_statistics.tsv",
                "output_dir": "./results/ml_out",
                "target_column": "class",
                "model": {"n_estimators": 25},
            }
        )
    )

    config = load_ml_config(str(config_path))

    assert config["target_column"] == "class"
    assert config["n_estimators"] == 25


def test_load_ml_config_defaults_min_samples_parameters(tmp_path):
    config_path = tmp_path / "ml_config_default_min_samples.json"
    config_path.write_text(
        json.dumps(
            {
                "input_path": "./results/summary_statistics.tsv",
                "output_dir": "./results/ml_out",
            }
        )
    )

    config = load_ml_config(str(config_path))

    assert config["min_samples_split"] == 2
    assert config["min_samples_leaf"] == 1


def test_load_ml_config_honors_overwrite_flag(tmp_path):
    config_path = tmp_path / "ml_config_overwrite_false.json"
    config_path.write_text(
        json.dumps(
            {
                "input_path": "./results/summary_statistics.tsv",
                "output_dir": "./results/ml_out",
                "overwrite": False,
            }
        )
    )

    config = load_ml_config(str(config_path))

    assert config["overwrite"] is False
