"""Normalize geNomad and ISEScan coordinates and compare predicted intervals."""

from __future__ import annotations

import re
from pathlib import Path

from io_utils import read_table

MOBILE_FIELDS = [
    "sample",
    "source",
    "kind",
    "element_id",
    "contig",
    "start",
    "end",
    "coordinate_status",
    "score",
    "fdr",
    "family",
]
ASSOCIATION_FIELDS = [
    "sample",
    "gene",
    "arg_contig",
    "arg_start",
    "arg_end",
    "source",
    "kind",
    "element_id",
    "same_contig",
    "relation",
    "distance_bp",
    "proven_membership",
]


def parse_genomad(path: Path, sample: str, kind: str, lengths: dict[str, int]) -> list[dict]:
    score = "plasmid_score" if kind == "plasmid" else "virus_score"
    rows = read_table(path, {"seq_name", "length", score})
    result = []
    for row in rows:
        name = row["seq_name"]
        # Only the documented terminal provirus suffix may be removed.
        match = re.fullmatch(r"(.+)\|provirus_(\d+)_(\d+)", name)
        contig = match[1] if match else name
        start, end = None, None
        state = "incompatible_identifier"
        if contig in lengths:
            state = "unavailable"
            if match:
                coordinates = row.get("coordinates", "")
                if re.fullmatch(r"\d+-\d+", coordinates):
                    start, end = map(int, coordinates.split("-"))
                    if (start, end) != (int(match[2]), int(match[3])):
                        state = "incompatible_coordinates"
                        start, end = None, None
            else:
                start, end = 1, lengths[contig]
            if start is not None:
                if 1 <= start <= end <= lengths[contig]:
                    state = "available"
                else:
                    start, end, state = None, None, "incompatible_coordinates"
        result.append(
            {
                "sample": sample,
                "source": "genomad",
                "kind": kind,
                "element_id": name,
                "contig": contig,
                "start": start,
                "end": end,
                "coordinate_status": state,
                "score": row[score],
                "fdr": row.get("fdr"),
                "family": None,
            }
        )
    return result


def parse_isescan(path: Path, sample: str, lengths: dict[str, int]) -> list[dict]:
    rows = read_table(path, {"seqID", "family", "isBegin", "isEnd"}, delimiter=",")
    result = []
    for i, row in enumerate(rows, 1):
        contig = row["seqID"]
        start, end = int(row["isBegin"]), int(row["isEnd"])
        state = "available"
        if contig not in lengths:
            state = "incompatible_identifier"
        elif not 1 <= start <= end <= lengths[contig]:
            state = "incompatible_coordinates"
        result.append(
            {
                "sample": sample,
                "source": "isescan",
                "kind": "IS",
                "element_id": f"IS_{i}",
                "contig": contig,
                "start": start if state == "available" else None,
                "end": end if state == "available" else None,
                "coordinate_status": state,
                "score": row.get("E-value"),
                "fdr": None,
                "family": row["family"],
            }
        )
    return result


def compare_intervals(a: tuple[int, int], b: tuple[int, int], threshold: int) -> tuple[str, int]:
    """Inputs are 1-based inclusive; distance counts bases strictly between intervals."""
    if threshold < 0 or any(start < 1 or end < start for start, end in (a, b)):
        raise ValueError("Invalid interval or distance threshold")
    if max(a[0], b[0]) <= min(a[1], b[1]):
        return "overlap", 0
    gap = max(a[0], b[0]) - min(a[1], b[1]) - 1
    return ("nearby" if gap <= threshold else "same_contig"), gap


def associate(
    args: list[dict],
    elements: list[dict],
    states: dict[str, str],
    lengths: dict[str, int],
    threshold: int,
) -> list[dict]:
    result = []
    for arg in args:
        if not arg["is_ARG"]:
            continue
        for source in ("genomad", "isescan"):
            base = {
                "sample": arg["sample"],
                "gene": arg["gene"],
                "arg_contig": arg["contig"],
                "arg_start": arg["start"],
                "arg_end": arg["end"],
                "source": source,
                "kind": None,
                "element_id": None,
                "same_contig": None,
                "distance_bp": None,
                "proven_membership": "not_established",
            }
            if states[source] != "success":
                result.append({**base, "relation": states[source]})
                continue
            if arg["contig"] not in lengths:
                result.append({**base, "relation": "incompatible_identifier"})
                continue
            if (
                arg["coordinate_status"] != "available"
                or not 1 <= arg["start"] <= arg["end"] <= lengths[arg["contig"]]
            ):
                result.append({**base, "relation": "arg_coordinates_unavailable"})
                continue
            candidates = [
                e for e in elements if e["source"] == source and e["contig"] == arg["contig"]
            ]
            if not candidates:
                incompatible = any(
                    e["source"] == source and e["coordinate_status"].startswith("incompatible")
                    for e in elements
                )
                relation = "incompatible_identifiers_present" if incompatible else "no_same_contig"
                result.append({**base, "same_contig": False, "relation": relation})
            for element in candidates:
                relation, gap = "coordinates_unavailable", None
                if element["coordinate_status"] == "available":
                    relation, gap = compare_intervals(
                        (arg["start"], arg["end"]), (element["start"], element["end"]), threshold
                    )
                result.append(
                    {
                        **base,
                        "kind": element["kind"],
                        "element_id": element["element_id"],
                        "same_contig": True,
                        "relation": relation,
                        "distance_bp": gap,
                    }
                )
    return result


def verify_isescan(work: Path, log: Path, sample: str, lengths: dict[str, int]) -> dict:
    """Require complete HMMER tables before accepting ISEScan findings or no-hits."""
    searches = [
        p
        for p in (work / "hmm").rglob("*")
        if p.is_file() and p.name.startswith(("clusters.faa.hmm.", "clusters.single.faa."))
    ]
    if len(searches) != 2 or any(not p.read_text().rstrip().endswith("# [ok]") for p in searches):
        raise ValueError("ISEScan HMMER searches are missing or incomplete")
    tables = list(work.rglob("assembly.fasta.csv"))
    if len(tables) > 1:
        raise ValueError("Ambiguous ISEScan native tables")
    if tables:
        parse_isescan(tables[0], sample, lengths)
        return {"native_table": str(tables[0].relative_to(work))}
    evidence = "No IS element was found for"
    if evidence not in log.read_text():
        raise ValueError("ISEScan produced no table and no verified no-hits message")
    return {"no_hits_evidence": evidence}
