rule preflight:
    input:
        manifest=config["manifest"],
        source=[str(p) for folder in ("scripts", "envs", "schemas")
                for p in (Path(workflow.basedir) / "workflow" / folder).glob("*") if p.is_file()],
        reads=[str(file) for record in SAMPLES.values() for file in
               (record.nanopore_reads, record.illumina_r1_reads, record.illumina_r2_reads) if file]
    output:
        directory(f"{OUT}/provenance")
    params:
        stage="preflight", outdir=lambda wc, output: str(output[0]), settings=config
    threads: allocation("preflight", "threads")
    resources:
        mem_mb=allocation("preflight", "mem_mb")
    log: "logs/preflight.log"
    benchmark: "benchmarks/preflight.tsv"
    conda: "../envs/python.yaml"
    script: "../scripts/dispatch.py"
