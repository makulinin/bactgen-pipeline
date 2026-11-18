import csv
import os

import pytest
from validate_manifest import COLUMNS, validate_manifest


def manifest(tmp_path, rows, columns=COLUMNS):
    for file in {v for row in rows for v in row[1:] if v and v.endswith((".fq", ".fastq.gz"))}:
        (tmp_path / file).write_bytes(b"@r\nACGT\n+\nIIII\n")
    path = tmp_path / "manifest.csv"
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(columns)
        writer.writerows(rows)
    return path


@pytest.mark.parametrize("ont", ["", "ont.fq"])
def test_valid(tmp_path, ont):
    path = manifest(tmp_path, [["S1", ont, "r1.fq", "r2.fq"]])
    sample = validate_manifest(path)["S1"]
    assert sample.illumina_r1_reads == (tmp_path / "r1.fq").resolve()
    assert bool(sample.nanopore_reads) == bool(ont)


@pytest.mark.parametrize("name", ["../S1", "bad id", "", "-bad", "S;rm", "a" * 101])
def test_invalid_id(tmp_path, name):
    with pytest.raises(ValueError, match="identifier"):
        validate_manifest(manifest(tmp_path, [[name, "", "r1.fq", "r2.fq"]]))


@pytest.mark.parametrize(
    "rows, message",
    [
        ([["S1", "", "r1.fq", "r2.fq"], ["S1", "", "r3.fq", "r4.fq"]], "duplicate"),
        ([["S1", "", "r1.fq", "r2.fq"], ["s1", "", "r3.fq", "r4.fq"]], "duplicate"),
        ([["S1", "ont.fq", "", ""]], "Nanopore-only"),
        ([["S1", "", "r1.fq", ""]], "R1 and R2"),
        ([["S1", "NA", "r1.fq", "r2.fq"]], "empty field"),
        ([["S1", "", "r1.fq", "r1.fq"]], "reused"),
        ([["S1", "", "r1.fq", "r2.fq"], ["S2", "", "r1.fq", "r3.fq"]], "reused"),
        ([["S1 ", "", "r1.fq", "r2.fq"]], "whitespace"),
        ([["S1", "", "r1.txt", "r2.fq"]], "extension"),
    ],
)
def test_rejected(tmp_path, rows, message):
    with pytest.raises(ValueError, match=message):
        validate_manifest(manifest(tmp_path, rows))


def test_columns(tmp_path):
    with pytest.raises(ValueError, match="columns"):
        validate_manifest(manifest(tmp_path, [], COLUMNS[:-1]))


def test_missing_file(tmp_path):
    path = manifest(tmp_path, [["S1", "", "r1.fq", "r2.fq"]])
    (tmp_path / "r1.fq").unlink()
    with pytest.raises(ValueError, match="missing"):
        validate_manifest(path)


def test_empty_file(tmp_path):
    path = manifest(tmp_path, [["S1", "", "r1.fq", "r2.fq"]])
    (tmp_path / "r2.fq").write_text("")
    with pytest.raises(ValueError, match="empty"):
        validate_manifest(path)


def test_hardlink_collision(tmp_path):
    path = manifest(tmp_path, [["S1", "", "r1.fq", "r2.fq"]])
    (tmp_path / "r2.fq").unlink()
    os.link(tmp_path / "r1.fq", tmp_path / "r2.fq")
    with pytest.raises(ValueError, match="reused"):
        validate_manifest(path)
