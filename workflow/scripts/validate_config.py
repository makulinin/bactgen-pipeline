"""Strict configuration validation, shared by CLI and Snakefile."""

from __future__ import annotations

import argparse
from pathlib import Path

import jsonschema
import yaml

STAGES = (
    "preflight",
    "fastqc",
    "fastp",
    "nanoplot",
    "multiqc",
    "unicycler",
    "quast",
    "checkm2",
    "bakta",
    "amrfinder",
    "genomad",
    "isescan",
    "summary",
)
MODULES = ("checkm2", "bakta", "amrfinder", "genomad", "isescan")


def obj(properties: dict) -> dict:
    return {
        "type": "object",
        "properties": properties,
        "required": list(properties),
        "additionalProperties": False,
    }


def number(low: float, high: float | None = None, integer: bool = False) -> dict:
    value = {"type": "integer" if integer else "number", "minimum": low}
    if high is not None:
        value["maximum"] = high
    return value


TEXT = {"type": "string"}
BOOL = {"type": "boolean"}
SCHEMA = obj(
    {
        "manifest": {"type": "string", "minLength": 1},
        "results_dir": {"type": "string", "minLength": 1},
        "modules": obj(dict.fromkeys(MODULES, BOOL)),
        "databases": obj({m: obj({"path": TEXT, "version": TEXT}) for m in MODULES[:-1]}),
        "resources": obj(
            {
                s: obj({"threads": number(1, integer=True), "mem_mb": number(256, integer=True)})
                for s in STAGES
            }
        ),
        "fastqc": obj({"enabled": BOOL}),
        "nanoplot": obj({"enabled": BOOL}),
        "fastp": obj(
            {
                "qualified_quality_phred": number(0, 93, True),
                "length_required": number(1, integer=True),
                "detect_adapter_for_pe": BOOL,
            }
        ),
        "unicycler": obj(
            {
                "mode": {"enum": ["normal", "bold", "conservative"]},
                "min_fasta_length": number(1, integer=True),
                "extra_args": {"type": "array", "items": TEXT},
            }
        ),
        "quast": obj({"min_contig": number(0, integer=True)}),
        "checkm2": obj(
            {
                "min_completeness": number(0, 100),
                "max_contamination": number(0),
                "exclude_failed": BOOL,
            }
        ),
        "bakta": obj(
            {
                "translation_table": {"enum": [4, 11, 25]},
                "genus": TEXT,
                "species": TEXT,
                "gram": {"enum": ["+", "-", "?"]},
                "min_contig_length": number(1, integer=True),
            }
        ),
        "amrfinder": obj({"plus": BOOL, "organism": TEXT}),
        "genomad": obj({"splits": number(1, integer=True), "cleanup": BOOL}),
        "isescan": obj({"remove_short_is": BOOL}),
        "coordinates": obj({"nearby_bp": number(0, integer=True)}),
    }
)


def validate_config(config: dict, databases: bool = False) -> None:
    jsonschema.validate(config, SCHEMA)
    output = Path(config["results_dir"])
    if output.is_absolute() or ".." in output.parts or str(output) in {".", ""}:
        raise ValueError("results_dir must be a dedicated relative subdirectory")
    if output.parts[0] in {"workflow", "tests", "resources", "data", "logs", "benchmarks"}:
        raise ValueError("results_dir overlaps a reserved source/input directory")
    if config["checkm2"]["exclude_failed"] and not config["modules"]["checkm2"]:
        raise ValueError("exclude_failed requires the CheckM2 module")
    # Additional assembly options are intentionally allowlisted to protect owned paths.
    args = config["unicycler"]["extra_args"]
    allowed = {
        "--min_kmer_frac",
        "--max_kmer_frac",
        "--kmer_count",
        "--min_anchor_seg_len",
        "--min_bridge_qual",
        "--linear_seqs",
        "--depth_filter",
    }
    if len(args) % 2 or any(args[i] not in allowed for i in range(0, len(args), 2)):
        raise ValueError("Unicycler extra_args must be allowlisted option/value pairs")
    for i in range(1, len(args), 2):
        if not args[i] or args[i].startswith("-"):
            raise ValueError("Unicycler additional values must be nonempty, nonnegative tokens")
    if databases:
        for module, db in config["databases"].items():
            if not config["modules"][module]:
                continue
            path = Path(db["path"])
            if not db["version"] or db["version"] == "SET_ME":
                raise ValueError(f"Record the exact {module} database release in config.yaml")
            if module == "checkm2":
                valid = path.is_file() and path.stat().st_size > 0
            else:
                valid = path.is_dir() and any(path.iterdir())
            if not valid:
                raise ValueError(f"Missing {module} database: {path}; see resources/README.md")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("config", type=Path)
    p.add_argument("--check-databases", action="store_true")
    a = p.parse_args()
    try:
        validate_config(yaml.safe_load(a.config.read_text()), a.check_databases)
    except (ValueError, OSError, jsonschema.ValidationError) as exc:
        p.exit(2, f"Configuration error: {exc}\n")
    print("Configuration valid")
