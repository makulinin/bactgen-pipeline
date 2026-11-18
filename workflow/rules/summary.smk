FUNCTIONAL_OUTPUTS = [f"{OUT}/{module}/{sample}" for module in
                      ("bakta", "amrfinder", "genomad", "isescan")
                      if config["modules"][module] for sample in IDS]

rule summarize:
    input:
        assemblies=expand(f"{OUT}/assembly/{{sample}}/assembly.fasta", sample=IDS),
        quast=expand(f"{OUT}/quast/{{sample}}", sample=IDS),
        checkm2=[f"{OUT}/checkm2"] if config["modules"]["checkm2"] else [],
        functional=FUNCTIONAL_OUTPUTS
    output:
        FINAL
    params:
        stage="summary", outdir=lambda wc, output: str(Path(output[0]).parent), samples=IDS, settings=config
    threads: allocation("summary", "threads")
    resources:
        mem_mb=allocation("summary", "mem_mb")
    log: "logs/summary.log"
    benchmark: "benchmarks/summary.tsv"
    conda: "../envs/python.yaml"
    script: "../scripts/dispatch.py"

rule annotation_all:
    input:
        FUNCTIONAL_OUTPUTS

rule summary_all:
    input:
        FINAL
