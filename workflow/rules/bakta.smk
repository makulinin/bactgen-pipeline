if config["modules"]["bakta"]:
    rule bakta:
        input:
            assembly=f"{OUT}/assembly/{{sample}}/assembly.fasta",
            quality=quality_input
        output:
            directory(f"{OUT}/bakta/{{sample}}")
        params:
            stage="bakta", outdir=lambda wc: f"{OUT}/bakta/{wc.sample}",
            settings=config["bakta"],
            quality_settings=config["checkm2"],
            database=config["databases"].get("bakta", {})
        threads: allocation("bakta", "threads")
        resources:
            mem_mb=allocation("bakta", "mem_mb")
        log: "logs/bakta/{sample}.log"
        benchmark: "benchmarks/bakta/{sample}.tsv"
        conda: "../envs/bakta.yaml"
        script: "../scripts/dispatch.py"
