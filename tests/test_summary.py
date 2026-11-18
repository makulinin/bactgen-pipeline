import csv
import json
import shutil

import pytest
from collect_annotations import parse_amrfinder, parse_bakta
from collect_assembly_stats import parse_checkm2, parse_quast, passes_quality
from conftest import DATA
from generate_summary import generate
from io_utils import fasta_lengths


def table(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_quast():
    row = parse_quast(DATA / "quast.tsv")
    assert row["genome_size"] == 2000
    assert row["N50"] == 1200
    assert row["L50"] == 1
    assert row["GC"] == 51.25


def test_checkm():
    row = parse_checkm2(DATA / "checkm2.tsv")["S1"]
    assert row["completeness"] == 98.2
    assert row["contamination"] == 1.1
    assert passes_quality(row, {"min_completeness": 90, "max_contamination": 5})
    assert not passes_quality(row, {"min_completeness": 99, "max_contamination": 5})


def test_annotations():
    assert parse_bakta(DATA / "annotation.gff3") == {"CDS": 1, "tRNA": 1, "rRNA": 1}
    rows = parse_amrfinder(DATA / "amrfinder.tsv", "S1")
    assert sum(row["is_ARG"] for row in rows) == 1
    assert parse_amrfinder(DATA / "amrfinder_empty.tsv", "S1") == []


def test_merge(result_tree):
    root, config = result_tree
    generate(config, ["S1"])
    row = table(root / "summary/samples_summary.tsv")[0]
    expected = {
        "assembly_status": "success",
        "contigs": "2",
        "genome_size": "2000",
        "completeness": "98.2",
        "CDS": "1",
        "ARG_count": "1",
        "AMR_point_mutations": "1",
        "plasmid_sequences": "1",
        "viral_sequences": "1",
        "IS_elements_count": "1",
    }
    assert all(row[k] == value for k, value in expected.items())
    assert "Completeness_Model_Used" in table(root / "summary/checkm2_normalized.tsv")[0]
    associations = table(root / "summary/arg_mobile_associations.tsv")
    assert {r["relation"] for r in associations} == {"overlap", "nearby"}
    assert all(r["proven_membership"] == "not_established" for r in associations)


@pytest.mark.parametrize("state", ["disabled", "missing", "skipped_quality"])
def test_missing_not_zero(result_tree, state):
    root, config = result_tree
    if state == "disabled":
        config["modules"]["amrfinder"] = False
    elif state == "missing":
        (root / "amrfinder/S1/status.json").unlink()
    else:
        (root / "amrfinder/S1/status.json").write_text(json.dumps({"status": state}))
    generate(config, ["S1"])
    row = table(root / "summary/samples_summary.tsv")[0]
    assert row["ARG_count"] == "NA"
    assert row["amrfinder_status"] == state


def test_real_empty_hits_zero(result_tree):
    root, config = result_tree
    shutil.copyfile(DATA / "amrfinder_empty.tsv", root / "amrfinder/S1/amrfinder.tsv")
    generate(config, ["S1"])
    assert table(root / "summary/resistance_summary.tsv")[0]["ARG_count"] == "0"


def test_corrupt_report_rejected(result_tree):
    root, config = result_tree
    (root / "quast/S1/report.tsv").write_text("wrong\theader\n")
    with pytest.raises(ValueError):
        generate(config, ["S1"])


@pytest.mark.parametrize("text", ["", ">c\n", "ACGT\n", ">c\nXYZ\n", ">c\nA\n>c\nA\n"])
def test_bad_fasta(tmp_path, text):
    path = tmp_path / "assembly.fasta"
    path.write_text(text)
    with pytest.raises(ValueError):
        fasta_lengths(path)


def test_gff_embedded_fasta(tmp_path):
    path = tmp_path / "annotation.gff3"
    path.write_text("##gff-version 3\n##FASTA\n>c1\nACGT\n")
    assert parse_bakta(path)["CDS"] == 0
