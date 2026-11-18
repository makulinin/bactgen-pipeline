"""Execute one Snakemake job using argument lists and promote validated outputs."""

from __future__ import annotations

import gzip
import hashlib
import json
import os
import shlex
import shutil
import subprocess
import sys
import tempfile
import time
import traceback
from pathlib import Path
from typing import Any


def execute(smk: Any) -> None:
    sys.path.insert(0, smk.scriptdir)
    from collect_annotations import parse_amrfinder, parse_bakta
    from collect_assembly_stats import parse_checkm2, parse_quast, passes_quality
    from collect_mobile_elements import parse_genomad
    from io_utils import fasta_lengths

    cfg, stage = smk.config, smk.params.stage
    target = Path(smk.params.outdir)
    log_path = Path(smk.log[0])
    log_path.parent.mkdir(parents=True, exist_ok=True)
    target.parent.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    sample = getattr(smk.wildcards, "sample", None)
    commands = []
    # Native environment records include build strings for reproducibility.
    metadata = []
    for file in sorted((Path(sys.prefix) / "conda-meta").glob("*.json")):
        record = json.loads(file.read_text())
        metadata.append({k: record.get(k) for k in ("name", "version", "build", "channel")})
    status = {
        "status": "success",
        "stage": stage,
        "sample": sample,
        "commands": commands,
        "environment": metadata,
    }

    def run(args: list, stdout: Path | None = None, cwd: Path | None = None) -> None:
        args = [str(value) for value in args]
        commands.append(args)
        log.write("$ " + shlex.join(args) + "\n")
        log.flush()
        env = dict(
            os.environ, OMP_NUM_THREADS=str(smk.threads), OPENBLAS_NUM_THREADS=str(smk.threads)
        )
        if stdout:
            with stdout.open("w") as handle:
                subprocess.run(args, stdout=handle, stderr=log, cwd=cwd, env=env, check=True)
        else:
            subprocess.run(args, stdout=log, stderr=log, cwd=cwd, env=env, check=True)

    def require(file: Path, nonempty: bool = True) -> None:
        if not file.is_file() or (nonempty and file.stat().st_size == 0):
            raise ValueError(f"Required output missing or empty: {file}")

    with (
        log_path.open("w") as log,
        tempfile.TemporaryDirectory(prefix=f".{target.name}.", dir=target.parent) as temporary,
    ):
        work = Path(temporary) / "output"
        try:
            if stage == "preflight":
                from validate_config import validate_config
                from validate_manifest import validate_manifest

                validate_config(cfg, databases=True)
                records = validate_manifest(cfg["manifest"])
                work.mkdir()
                shutil.copy2(cfg["manifest"], work / "manifest.csv")
                (work / "config.json").write_text(json.dumps(cfg, indent=2))
                inputs = {
                    name: {key: str(value) if value else None for key, value in vars(rec).items()}
                    for name, rec in records.items()
                }
                (work / "resolved_manifest.json").write_text(json.dumps(inputs, indent=2))
                status["environment_specs_sha256"] = {
                    file.name: hashlib.sha256(file.read_bytes()).hexdigest()
                    for file in Path(smk.scriptdir).parent.joinpath("envs").glob("*.yaml")
                }
            elif stage == "summary":
                from generate_summary import generate

                # Write each table atomically; Snakemake removes outputs on job failure.
                generate(cfg, list(smk.params.samples))
                log.write("Summary tables generated\n")
                return
            else:
                excluded = False
                if (
                    stage in {"bakta", "amrfinder", "genomad", "isescan"}
                    and cfg["checkm2"]["exclude_failed"]
                ):
                    report = Path(smk.input.quality) / "quality_report.tsv"
                    row = parse_checkm2(report)[sample]
                    excluded = not passes_quality(row, cfg["checkm2"])
                if excluded:
                    work.mkdir()
                    status["status"] = "skipped_quality"
                elif stage == "fastp":
                    work.mkdir()
                    settings = cfg["fastp"]
                    args = [
                        "fastp",
                        "--in1",
                        smk.input.r1,
                        "--in2",
                        smk.input.r2,
                        "--out1",
                        work / "R1.fastq.gz",
                        "--out2",
                        work / "R2.fastq.gz",
                        "--json",
                        work / "fastp.json",
                        "--html",
                        work / "fastp.html",
                        "--thread",
                        smk.threads,
                        "--qualified_quality_phred",
                        settings["qualified_quality_phred"],
                        "--length_required",
                        settings["length_required"],
                    ]
                    if settings["detect_adapter_for_pe"]:
                        args.append("--detect_adapter_for_pe")
                    run(args)
                    for file in ("R1.fastq.gz", "R2.fastq.gz", "fastp.json", "fastp.html"):
                        require(work / file)
                    json.loads((work / "fastp.json").read_text())
                    # Ensure the cleaned files are readable and contain at least one full record.
                    for mate in ("R1", "R2"):
                        with gzip.open(work / f"{mate}.fastq.gz", "rt") as handle:
                            record = [handle.readline().rstrip("\n") for _ in range(4)]
                        if (
                            not record[0].startswith("@")
                            or not record[2].startswith("+")
                            or not record[1]
                            or len(record[1]) != len(record[3])
                        ):
                            raise ValueError("fastp produced empty or invalid cleaned FASTQ")
                elif stage == "fastqc":
                    work.mkdir()
                    source = Path(smk.input.reads)
                    suffix = ".fastq.gz" if str(source).lower().endswith(".gz") else ".fastq"
                    prefix = f"{sample}_{target.name}"
                    link = Path(temporary) / (prefix + suffix)
                    link.symlink_to(source.resolve())
                    run(["fastqc", "--threads", smk.threads, "--outdir", work, link])
                    require(work / f"{prefix}_fastqc.html")
                    require(work / f"{prefix}_fastqc.zip")
                elif stage == "nanoplot":
                    run(
                        [
                            "NanoPlot",
                            "--fastq",
                            smk.input.reads,
                            "--threads",
                            smk.threads,
                            "--outdir",
                            work,
                        ]
                    )
                    require(work / "NanoPlot-report.html")
                elif stage == "multiqc":
                    run(
                        [
                            "multiqc",
                            "--dirs",
                            "--dirs-depth",
                            "2",
                            "--force",
                            "--outdir",
                            work,
                            "--filename",
                            "multiqc_report.html",
                            *smk.input.reports,
                        ]
                    )
                    require(work / "multiqc_report.html")
                elif stage == "unicycler":
                    settings = cfg["unicycler"]
                    args = [
                        "unicycler",
                        "-1",
                        smk.input.r1,
                        "-2",
                        smk.input.r2,
                        "-o",
                        work,
                        "--threads",
                        smk.threads,
                        "--mode",
                        settings["mode"],
                        "--min_fasta_length",
                        settings["min_fasta_length"],
                    ]
                    if smk.params.nanopore:
                        args.extend(["-l", smk.params.nanopore])
                    args.extend(settings["extra_args"])
                    run(args)
                    fasta_lengths(work / "assembly.fasta")
                    require(work / "unicycler.log")
                    status["gfa_available"] = (work / "assembly.gfa").is_file()
                elif stage == "quast":
                    run(
                        [
                            "quast.py",
                            smk.input.assembly,
                            "--output-dir",
                            work,
                            "--threads",
                            smk.threads,
                            "--min-contig",
                            cfg["quast"]["min_contig"],
                        ]
                    )
                    parse_quast(work / "report.tsv")
                elif stage == "checkm2":
                    bins = Path(temporary) / "bins"
                    bins.mkdir()
                    for name, source in zip(smk.params.samples, smk.input.assemblies):
                        (bins / f"{name}.fasta").symlink_to(Path(source).resolve())
                    run(
                        [
                            "checkm2",
                            "predict",
                            "--input",
                            bins,
                            "--output-directory",
                            work,
                            "--threads",
                            smk.threads,
                            "--extension",
                            "fasta",
                            "--database_path",
                            cfg["databases"]["checkm2"]["path"],
                        ]
                    )
                    parsed = parse_checkm2(work / "quality_report.tsv")
                    if set(parsed) != set(smk.params.samples):
                        raise ValueError("CheckM2 output sample set differs from input")
                elif stage == "bakta":
                    settings = cfg["bakta"]
                    args = [
                        "bakta",
                        "--db",
                        cfg["databases"]["bakta"]["path"],
                        "--output",
                        work,
                        "--prefix",
                        "annotation",
                        "--threads",
                        smk.threads,
                        "--keep-contig-headers",
                        "--translation-table",
                        settings["translation_table"],
                        "--min-contig-length",
                        settings["min_contig_length"],
                        "--gram",
                        settings["gram"],
                    ]
                    for key in ("genus", "species"):
                        if settings[key]:
                            args.extend(["--" + key, settings[key]])
                    run([*args, smk.input.assembly])
                    for ext in ("gff3", "gbff", "faa", "ffn", "tsv", "json"):
                        require(work / f"annotation.{ext}", nonempty=ext not in {"faa", "ffn"})
                    parse_bakta(work / "annotation.gff3")
                    json.loads((work / "annotation.json").read_text())
                elif stage == "amrfinder":
                    work.mkdir()
                    args = [
                        "amrfinder",
                        "-n",
                        smk.input.assembly,
                        "--database",
                        cfg["databases"]["amrfinder"]["path"],
                        "--threads",
                        smk.threads,
                        "--output",
                        work / "amrfinder.tsv",
                    ]
                    if cfg["amrfinder"]["plus"]:
                        args.append("--plus")
                    if cfg["amrfinder"]["organism"]:
                        args.extend(["--organism", cfg["amrfinder"]["organism"]])
                    run(args)
                    parse_amrfinder(work / "amrfinder.tsv", sample)
                elif stage == "genomad":
                    args = [
                        "genomad",
                        "end-to-end",
                        "--threads",
                        smk.threads,
                        "--splits",
                        cfg["genomad"]["splits"],
                    ]
                    if cfg["genomad"]["cleanup"]:
                        args.append("--cleanup")
                    run([*args, smk.input.assembly, work, cfg["databases"]["genomad"]["path"]])
                    for kind in ("plasmid", "virus"):
                        parse_genomad(
                            work / "assembly_summary" / f"assembly_{kind}_summary.tsv",
                            sample,
                            kind,
                            fasta_lengths(Path(smk.input.assembly)),
                        )
                elif stage == "isescan":
                    args = [
                        "isescan.py",
                        "--seqfile",
                        Path(smk.input.assembly).resolve(),
                        "--output",
                        work.resolve(),
                        "--nthread",
                        smk.threads,
                    ]
                    if cfg["isescan"]["remove_short_is"]:
                        args.append("--removeShortIS")
                    run(args, cwd=Path(temporary))
                    from collect_mobile_elements import verify_isescan

                    status.update(
                        verify_isescan(
                            work, log_path, sample, fasta_lengths(Path(smk.input.assembly))
                        )
                    )
                else:
                    raise ValueError(f"Unknown stage: {stage}")
            status["elapsed_seconds"] = round(time.monotonic() - started, 3)
            (work / "status.json").write_text(json.dumps(status, indent=2) + "\n")
            if target.exists():
                shutil.rmtree(target)
            os.replace(work, target)
        except Exception:
            log.write(traceback.format_exc())
            log.flush()
            raise
