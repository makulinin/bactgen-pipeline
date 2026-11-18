if config["modules"]["checkm2"]:
    rule checkm2:
        input:
            assemblies=expand(f"{OUT}/assembly/{{sample}}/assembly.fasta", sample=IDS)
        output:
            directory(f"{OUT}/checkm2")
        params:
            stage="checkm2", outdir=lambda wc, output: str(output[0]), samples=IDS,
            settings=config["checkm2"], database=config["databases"]["checkm2"]
        threads: allocation("checkm2", "threads")
        resources:
            mem_mb=allocation("checkm2", "mem_mb")
        log: "logs/checkm2.log"
        benchmark: "benchmarks/checkm2.tsv"
        conda: "../envs/checkm2.yaml"
        script: "../scripts/dispatch.py"
