"""Parse Bakta GFF3 and AMRFinderPlus 3.12 native tabular results."""

from __future__ import annotations

from collections import Counter
from pathlib import Path

from io_utils import read_table

AMR_FIELDS = [
    "sample",
    "contig",
    "start",
    "end",
    "strand",
    "gene",
    "element_type",
    "element_subtype",
    "class",
    "subclass",
    "method",
    "is_ARG",
    "coordinate_status",
]


def parse_bakta(path: Path) -> dict:
    counts: Counter = Counter()
    with path.open() as handle:
        for line in handle:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                raise ValueError(f"Malformed GFF3 line in {path}")
            counts[parts[2]] += 1
    return {"CDS": counts["CDS"], "tRNA": counts["tRNA"], "rRNA": counts["rRNA"]}


def parse_amrfinder(path: Path, sample: str) -> list[dict]:
    required = {
        "Contig id",
        "Start",
        "Stop",
        "Strand",
        "Gene symbol",
        "Element type",
        "Element subtype",
        "Class",
        "Subclass",
        "Method",
    }
    rows = read_table(path, required)
    result = []
    for row in rows:
        start, end = row["Start"], row["Stop"]
        valid = start.isdigit() and end.isdigit() and int(start) > 0 and int(end) > 0
        bounds = sorted([int(start), int(end)]) if valid else [None, None]
        result.append(
            {
                "sample": sample,
                "contig": row["Contig id"],
                "start": bounds[0],
                "end": bounds[1],
                "strand": row["Strand"],
                "gene": row["Gene symbol"],
                "element_type": row["Element type"],
                "element_subtype": row["Element subtype"],
                "class": row["Class"],
                "subclass": row["Subclass"],
                "method": row["Method"],
                "is_ARG": row["Element type"] == "AMR" and row["Element subtype"] == "AMR",
                "coordinate_status": "available" if valid else "unavailable",
            }
        )
    return result
