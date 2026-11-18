import pytest
from jsonschema import ValidationError
from validate_config import validate_config


def test_valid(config):
    validate_config(config)


@pytest.mark.parametrize("change", ["threads", "unknown", "mode", "gate", "extra", "out"])
def test_invalid(config, change):
    if change == "threads":
        config["resources"]["fastp"]["threads"] = 0
    elif change == "unknown":
        config["magic"] = True
    elif change == "mode":
        config["unicycler"]["mode"] = "made-up"
    elif change == "gate":
        config["modules"]["checkm2"] = False
        config["checkm2"]["exclude_failed"] = True
    elif change == "extra":
        config["unicycler"]["extra_args"] = ["--out", "somewhere"]
    else:
        config["results_dir"] = "."
    with pytest.raises((ValueError, ValidationError)):
        validate_config(config)


def test_database_preflight(config):
    with pytest.raises(ValueError, match="release"):
        validate_config(config, databases=True)
