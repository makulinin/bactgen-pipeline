rule assemble:
    input:
        r1=f"{OUT}/qc/fastp/{{sample}}/R1.fastq.gz",
        r2=f"{OUT}/qc/fastp/{{sample}}/R2.fastq.gz",
        nanopore=lambda wc: [str(SAMPLES[wc.sample].nanopore_reads)] if
                            SAMPLES[wc.sample].nanopore_reads else []
    output:
        fasta=ensure(f"{OUT}/assembly/{{sample}}/assembly.fasta", non_empty=True),
        status=f"{OUT}/assembly/{{sample}}/status.json"
    params:
        stage="unicycler", outdir=lambda wc: f"{OUT}/assembly/{wc.sample}",
        nanopore=lambda wc: str(SAMPLES[wc.sample].nanopore_reads) if
                             SAMPLES[wc.sample].nanopore_reads else "",
        settings=config["unicycler"]
    threads: allocation("unicycler", "threads")
    resources:
        mem_mb=allocation("unicycler", "mem_mb")
    log: "logs/unicycler/{sample}.log"
    benchmark: "benchmarks/unicycler/{sample}.tsv"
    conda: "../envs/unicycler.yaml"
    script: "../scripts/dispatch.py"

rule assembly_all:
    input:
        expand(f"{OUT}/assembly/{{sample}}/assembly.fasta", sample=IDS)
