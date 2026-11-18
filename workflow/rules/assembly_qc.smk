rule quast:
    input:
        assembly=f"{OUT}/assembly/{{sample}}/assembly.fasta"
    output:
        directory(f"{OUT}/quast/{{sample}}")
    params:
        stage="quast", outdir=lambda wc: f"{OUT}/quast/{wc.sample}", settings=config["quast"]
    threads: allocation("quast", "threads")
    resources:
        mem_mb=allocation("quast", "mem_mb")
    log: "logs/quast/{sample}.log"
    benchmark: "benchmarks/quast/{sample}.tsv"
    conda: "../envs/quast.yaml"
    script: "../scripts/dispatch.py"

rule assembly_qc_all:
    input:
        expand(f"{OUT}/quast/{{sample}}", sample=IDS),
        [f"{OUT}/checkm2"] if config["modules"]["checkm2"] else []
