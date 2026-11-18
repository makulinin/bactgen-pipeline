"""Collect completed analyses; status fields distinguish missing, skipped and zero."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import yaml
from collect_annotations import AMR_FIELDS, parse_amrfinder, parse_bakta
from collect_assembly_stats import parse_checkm2, parse_quast, passes_quality
from collect_mobile_elements import (
    ASSOCIATION_FIELDS,
    MOBILE_FIELDS,
    associate,
    parse_genomad,
    parse_isescan,
)
from io_utils import analysis_status, fasta_lengths, write_tsv
from validate_manifest import validate_manifest

ASSEMBLY_FIELDS = [
    "sample",
    "assembly_status",
    "quast_status",
    "genome_size",
    "contigs",
    "N50",
    "L50",
    "GC",
    "max_contig_length",
    "checkm2_status",
    "completeness",
    "contamination",
    "quality_pass",
]
ANNOTATION_FIELDS = ["sample", "bakta_status", "CDS", "tRNA", "rRNA"]
RESISTANCE_FIELDS = [
    "sample",
    "amrfinder_status",
    "ARG_count",
    "AMR_point_mutations",
    "AMRFinder_total_hits",
]
MOBILE_SUMMARY_FIELDS = [
    "sample",
    "genomad_status",
    "isescan_status",
    "plasmid_sequences",
    "viral_sequences",
    "IS_elements_count",
]
SUMMARY_NAMES = [
    "assembly_summary",
    "annotation_summary",
    "resistance_summary",
    "mobile_elements_summary",
    "samples_summary",
]


def generate(config: dict, samples: list[str]) -> None:
    root = Path(config["results_dir"])
    out = root / "summary"
    module = config["modules"]
    checkm_path = root / "checkm2" / "quality_report.tsv"
    checkm = parse_checkm2(checkm_path) if module["checkm2"] and checkm_path.exists() else {}
    if checkm and set(checkm) != set(samples):
        raise ValueError("CheckM2 report and manifest sample sets disagree")
    rows, all_args, all_mobile, associations, normalized_checkm = [], [], [], [], []
    for sample in samples:
        row = {"sample": sample}
        fasta = root / "assembly" / sample / "assembly.fasta"
        lengths = fasta_lengths(fasta) if fasta.exists() else {}
        row["assembly_status"] = "success" if lengths else "missing"
        quast = root / "quast" / sample / "report.tsv"
        row["quast_status"] = "success" if quast.exists() else "missing"
        if quast.exists():
            metrics = parse_quast(quast)
            metrics.pop("quast_metrics")
            row.update(metrics)
        row["checkm2_status"] = "disabled" if not module["checkm2"] else "missing"
        if sample in checkm:
            metrics = {k: v for k, v in checkm[sample].items() if k != "checkm2_metrics"}
            row.update(metrics)
            row["checkm2_status"] = "success"
            row["quality_pass"] = passes_quality(checkm[sample], config["checkm2"])
            normalized_checkm.append(
                {
                    **checkm[sample]["checkm2_metrics"],
                    "sample": sample,
                    **metrics,
                    "quality_pass": row["quality_pass"],
                }
            )
        states = {
            m: analysis_status(root, m, sample, module[m])
            for m in ("bakta", "amrfinder", "genomad", "isescan")
        }
        row.update({f"{m}_status": state for m, state in states.items()})
        if states["bakta"] == "success":
            row.update(parse_bakta(root / "bakta" / sample / "annotation.gff3"))
        args, elements = [], []
        if states["amrfinder"] == "success":
            args = parse_amrfinder(root / "amrfinder" / sample / "amrfinder.tsv", sample)
            row["ARG_count"] = sum(a["is_ARG"] for a in args)
            row["AMR_point_mutations"] = sum(
                a["element_type"] == "AMR" and a["element_subtype"] == "POINT" for a in args
            )
            row["AMRFinder_total_hits"] = len(args)
        if states["genomad"] == "success":
            for kind, field in (("plasmid", "plasmid_sequences"), ("virus", "viral_sequences")):
                path = (
                    root / "genomad" / sample / "assembly_summary" / f"assembly_{kind}_summary.tsv"
                )
                found = parse_genomad(path, sample, kind, lengths)
                elements.extend(found)
                row[field] = len(found)
        if states["isescan"] == "success":
            status = json.loads((root / "isescan" / sample / "status.json").read_text())
            if status.get("native_table"):
                found = parse_isescan(
                    root / "isescan" / sample / status["native_table"], sample, lengths
                )
            elif status.get("no_hits_evidence"):
                found = []
            else:
                raise ValueError(f"ISEScan missing table and no-hits evidence for {sample}")
            row["IS_elements_count"] = len(found)
            elements.extend(found)
        all_args.extend(args)
        all_mobile.extend(elements)
        associations.extend(
            associate(args, elements, states, lengths, config["coordinates"]["nearby_bp"])
        )
        rows.append(row)
    fields = (
        ASSEMBLY_FIELDS + ANNOTATION_FIELDS[1:] + RESISTANCE_FIELDS[1:] + MOBILE_SUMMARY_FIELDS[1:]
    )
    for name, columns in zip(
        SUMMARY_NAMES,
        [ASSEMBLY_FIELDS, ANNOTATION_FIELDS, RESISTANCE_FIELDS, MOBILE_SUMMARY_FIELDS, fields],
    ):
        write_tsv(out / f"{name}.tsv", rows, columns)
    write_tsv(out / "resistance_genes.tsv", all_args, AMR_FIELDS)
    write_tsv(out / "mobile_elements.tsv", all_mobile, MOBILE_FIELDS)
    write_tsv(out / "arg_mobile_associations.tsv", associations, ASSOCIATION_FIELDS)
    write_tsv(
        out / "checkm2_normalized.tsv",
        normalized_checkm,
        ["sample", "completeness", "contamination", "quality_pass"]
        + sorted(
            {key for row in normalized_checkm for key in row}
            - {
                "sample",
                "completeness",
                "contamination",
                "quality_pass",
                "Name",
                "Completeness",
                "Contamination",
            }
        ),
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--config", type=Path, default=Path("config.yaml"))
    a = p.parse_args()
    cfg = yaml.safe_load(a.config.read_text())
    generate(cfg, list(validate_manifest(cfg["manifest"])))
