if config["modules"]["amrfinder"]:
    rule amrfinder:
        input:
            assembly=f"{OUT}/assembly/{{sample}}/assembly.fasta",
            quality=quality_input
        output:
            directory(f"{OUT}/amrfinder/{{sample}}")
        params:
            stage="amrfinder", outdir=lambda wc: f"{OUT}/amrfinder/{wc.sample}",
            settings=config["amrfinder"],
            quality_settings=config["checkm2"],
            database=config["databases"].get("amrfinder", {})
        threads: allocation("amrfinder", "threads")
        resources:
            mem_mb=allocation("amrfinder", "mem_mb")
        log: "logs/amrfinder/{sample}.log"
        benchmark: "benchmarks/amrfinder/{sample}.tsv"
        conda: "../envs/amrfinder.yaml"
        script: "../scripts/dispatch.py"
