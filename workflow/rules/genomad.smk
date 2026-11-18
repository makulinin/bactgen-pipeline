if config["modules"]["genomad"]:
    rule genomad:
        input:
            assembly=f"{OUT}/assembly/{{sample}}/assembly.fasta",
            quality=quality_input
        output:
            directory(f"{OUT}/genomad/{{sample}}")
        params:
            stage="genomad", outdir=lambda wc: f"{OUT}/genomad/{wc.sample}",
            settings=config["genomad"],
            quality_settings=config["checkm2"],
            database=config["databases"].get("genomad", {})
        threads: allocation("genomad", "threads")
        resources:
            mem_mb=allocation("genomad", "mem_mb")
        log: "logs/genomad/{sample}.log"
        benchmark: "benchmarks/genomad/{sample}.tsv"
        conda: "../envs/genomad.yaml"
        script: "../scripts/dispatch.py"
