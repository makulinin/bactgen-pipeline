import pytest
from collect_mobile_elements import associate, compare_intervals, parse_genomad, parse_isescan
from conftest import DATA


@pytest.mark.parametrize(
    "a,b,limit,expected",
    [
        ((10, 20), (20, 30), 0, ("overlap", 0)),
        ((10, 20), (21, 30), 0, ("nearby", 0)),
        ((10, 20), (22, 30), 1, ("nearby", 1)),
        ((10, 20), (22, 30), 0, ("same_contig", 1)),
        ((30, 40), (10, 20), 9, ("nearby", 9)),
    ],
)
def test_intervals(a, b, limit, expected):
    assert compare_intervals(a, b, limit) == expected
    assert compare_intervals(b, a, limit) == expected


def test_provirus_coordinates():
    row = parse_genomad(DATA / "virus.tsv", "S1", "virus", {"c1": 1200})[0]
    assert (row["contig"], row["start"], row["end"]) == ("c1", 50, 200)
    assert row["coordinate_status"] == "available"
    assert (
        parse_genomad(DATA / "virus.tsv", "S1", "virus", {"c2": 1200})[0]["coordinate_status"]
        == "incompatible_identifier"
    )


def test_is_coordinates():
    row = parse_isescan(DATA / "isescan.csv", "S1", {"c1": 1200})[0]
    assert (row["start"], row["end"], row["family"]) == (103, 190, "IS3")


def test_unknown_arg_coordinates():
    args = [
        {
            "sample": "S1",
            "gene": "arg",
            "contig": "c1",
            "start": None,
            "end": None,
            "is_ARG": True,
            "coordinate_status": "unavailable",
        }
    ]
    result = associate(args, [], {"genomad": "success", "isescan": "disabled"}, {"c1": 100}, 10)
    assert [r["relation"] for r in result] == ["arg_coordinates_unavailable", "disabled"]


def test_invalid_intervals():
    with pytest.raises(ValueError):
        compare_intervals((5, 1), (1, 2), 10)


@pytest.mark.parametrize("complete,empty_message", [(True, True), (False, True), (True, False)])
def test_isescan_no_hits_evidence(tmp_path, complete, empty_message):
    from collect_mobile_elements import verify_isescan

    hmm = tmp_path / "hmm"
    hmm.mkdir()
    for prefix in ("clusters.faa.hmm", "clusters.single.faa"):
        (hmm / (prefix + ".assembly.fasta.faa")).write_text(
            "# HMMER fixture\n" + ("# [ok]\n" if complete else "")
        )
    log = tmp_path / "isescan.log"
    log.write_text("No IS element was found for fixture" if empty_message else "failed")
    if complete and empty_message:
        assert verify_isescan(tmp_path, log, "S1", {})["no_hits_evidence"]
    else:
        with pytest.raises(ValueError):
            verify_isescan(tmp_path, log, "S1", {})
