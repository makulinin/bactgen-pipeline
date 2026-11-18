"""Validate manifests; resolve reads relative to the manifest, never the shell."""

from __future__ import annotations

import argparse
import csv
import re
from dataclasses import dataclass
from pathlib import Path

COLUMNS = ["sample", "nanopore_reads", "illumina_r1_reads", "illumina_r2_reads"]
EXTENSIONS = (".fastq", ".fq", ".fastq.gz", ".fq.gz")


@dataclass(frozen=True)
class Sample:
    sample: str
    nanopore_reads: Path | None
    illumina_r1_reads: Path
    illumina_r2_reads: Path


def validate_manifest(path: str | Path) -> dict[str, Sample]:
    path = Path(path).resolve()
    samples: dict[str, Sample] = {}
    used: dict[tuple[int, int], str] = {}
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != COLUMNS:
            raise ValueError(f"Manifest columns must be exactly {','.join(COLUMNS)} in this order")
        for line, row in enumerate(reader, 2):
            if None in row or any(v is None for v in row.values()):
                raise ValueError(f"Line {line}: expected exactly four fields")
            if any(v != v.strip() for v in row.values()):
                raise ValueError(f"Line {line}: surrounding whitespace is not allowed")
            name = row["sample"]
            if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", name) or len(name) > 100:
                raise ValueError(f"Line {line}: invalid sample identifier {name!r}")
            if name.casefold() in {s.casefold() for s in samples}:
                raise ValueError(f"Line {line}: duplicate sample (case insensitive): {name}")
            if not row["illumina_r1_reads"] or not row["illumina_r2_reads"]:
                raise ValueError(
                    f"{name}: paired Illumina R1 and R2 are required; Nanopore-only unsupported"
                )
            reads: list[Path | None] = []
            for column in COLUMNS[1:]:
                value = row[column]
                if value.lower() in {"na", "nan", "null", "none", "-"}:
                    raise ValueError(f"{name}: use an empty field, not {value!r}")
                if not value:
                    reads.append(None)
                    continue
                file = (path.parent / value).resolve()
                if not str(file).lower().endswith(EXTENSIONS):
                    raise ValueError(f"{name}: unsupported FASTQ extension: {file}")
                if not file.is_file() or not file.stat().st_size:
                    raise ValueError(f"{name}: missing or empty reads: {file}")
                stat = file.stat()
                identity = (stat.st_dev, stat.st_ino)
                if identity in used:
                    raise ValueError(f"{name}: reads reused by {used[identity]}: {file}")
                used[identity] = name
                reads.append(file)
            samples[name] = Sample(name, *reads)
    if not samples:
        raise ValueError("Manifest contains no samples; populate manifest.csv before running")
    return samples


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    args = parser.parse_args()
    try:
        result = validate_manifest(args.manifest)
    except (ValueError, OSError) as exc:
        parser.exit(2, f"Manifest error: {exc}\n")
    print(f"Validated {len(result)} samples")
