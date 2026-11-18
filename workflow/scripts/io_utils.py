"""Small strict readers and atomic writers used by all result collectors."""

from __future__ import annotations

import csv
import json
import math
import os
import re
import tempfile
from pathlib import Path
from typing import Any, Iterable


def write_tsv(path: Path, rows: Iterable[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(fd, "w", newline="") as handle:
            writer = csv.DictWriter(handle, fields, delimiter="\t", extrasaction="ignore")
            writer.writeheader()
            for row in rows:
                writer.writerow({k: "NA" if row.get(k) is None else row.get(k) for k in fields})
        os.replace(name, path)
    finally:
        Path(name).unlink(missing_ok=True)


def read_table(path: Path, required: set[str], delimiter: str = "\t") -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f"{path}: missing columns {required - set(reader.fieldnames or [])}")
        rows = list(reader)
        if any(None in row or any(v is None for v in row.values()) for row in rows):
            raise ValueError(f"{path}: malformed table")
        return rows


def numeric(value: Any) -> float | None:
    if value is None or str(value).strip() in {"", "NA", "N/A", "-"}:
        return None
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"Nonfinite numeric result: {value}")
    return result


def fasta_lengths(path: Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    name: str | None = None
    with path.open() as handle:
        for line in handle:
            text = line.strip()
            if not text:
                continue
            if text.startswith(">"):
                if name is not None and lengths[name] == 0:
                    raise ValueError(f"{path}: empty FASTA record {name}")
                parts = text[1:].split()
                if not parts or parts[0] in lengths:
                    raise ValueError(f"{path}: empty or duplicate FASTA identifier")
                name = parts[0]
                lengths[name] = 0
            elif name is None or not re.fullmatch(r"[ACGTRYSWKMBDHVNacgtryswkmbdhvn]+", text):
                raise ValueError(f"{path}: invalid nucleotide FASTA")
            else:
                lengths[name] += len(text)
    if not lengths or lengths[name] == 0:
        raise ValueError(f"{path}: empty nucleotide FASTA")
    return lengths


def analysis_status(root: Path, module: str, sample: str, enabled: bool) -> str:
    if not enabled:
        return "disabled"
    file = root / module / sample / "status.json"
    if not file.exists():
        return "missing"
    value = json.loads(file.read_text())["status"]
    if value not in {"success", "skipped_quality"}:
        raise ValueError(f"Unexpected analysis status: {value}")
    return value
