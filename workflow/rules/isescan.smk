if config["modules"]["isescan"]:
    rule isescan:
        input:
            assembly=f"{OUT}/assembly/{{sample}}/assembly.fasta",
            quality=quality_input
        output:
            directory(f"{OUT}/isescan/{{sample}}")
        params:
            stage="isescan", outdir=lambda wc: f"{OUT}/isescan/{wc.sample}",
            settings=config["isescan"],
            quality_settings=config["checkm2"],
            database=config["databases"].get("isescan", {})
        threads: allocation("isescan", "threads")
        resources:
            mem_mb=allocation("isescan", "mem_mb")
        log: "logs/isescan/{sample}.log"
        benchmark: "benchmarks/isescan/{sample}.tsv"
        conda: "../envs/isescan.yaml"
        script: "../scripts/dispatch.py"
