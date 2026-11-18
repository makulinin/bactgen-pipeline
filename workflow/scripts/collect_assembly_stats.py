"""Parse native QUAST and CheckM2 reports without replacing unknowns with zero."""

from __future__ import annotations

import csv
from pathlib import Path

from io_utils import numeric, read_table

QUAST = {
    "Total length": "genome_size",
    "# contigs": "contigs",
    "N50": "N50",
    "L50": "L50",
    "GC (%)": "GC",
    "Largest contig": "max_contig_length",
}


def parse_quast(path: Path) -> dict:
    with path.open(newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    if not rows or rows[0][0] != "Assembly" or any(len(row) != 2 for row in rows):
        raise ValueError(f"{path}: expected a single-assembly QUAST report")
    raw = {row[0]: row[1] for row in rows[1:]}
    if not set(QUAST).issubset(raw):
        raise ValueError(f"{path}: incomplete QUAST metrics")
    result = {dest: numeric(raw[source]) for source, dest in QUAST.items()}
    for field in set(QUAST.values()) - {"GC"}:
        if result[field] is not None:
            if not result[field].is_integer():
                raise ValueError(f"{path}: expected integer metric {field}")
            result[field] = int(result[field])
    result["quast_metrics"] = raw
    return result


def parse_checkm2(path: Path) -> dict[str, dict]:
    rows = read_table(path, {"Name", "Completeness", "Contamination"})
    result = {}
    for row in rows:
        name = row["Name"]
        if name in result:
            raise ValueError(f"{path}: duplicate CheckM2 sample {name}")
        completeness = numeric(row["Completeness"])
        contamination = numeric(row["Contamination"])
        if (
            completeness is None
            or contamination is None
            or not 0 <= completeness <= 100
            or contamination < 0
        ):
            raise ValueError(f"{path}: invalid CheckM2 predictions for {name}")
        result[name] = {
            "completeness": completeness,
            "contamination": contamination,
            "checkm2_metrics": row,
        }
    return result


def passes_quality(row: dict, settings: dict) -> bool:
    return (
        row["completeness"] >= settings["min_completeness"]
        and row["contamination"] <= settings["max_contamination"]
    )
