rule fastp:
    input:
        preflight=f"{OUT}/provenance",
        r1=lambda wc: str(SAMPLES[wc.sample].illumina_r1_reads),
        r2=lambda wc: str(SAMPLES[wc.sample].illumina_r2_reads)
    output:
        r1=f"{OUT}/qc/fastp/{{sample}}/R1.fastq.gz",
        r2=f"{OUT}/qc/fastp/{{sample}}/R2.fastq.gz",
        json=f"{OUT}/qc/fastp/{{sample}}/fastp.json",
        html=f"{OUT}/qc/fastp/{{sample}}/fastp.html"
    params:
        stage="fastp", outdir=lambda wc: f"{OUT}/qc/fastp/{wc.sample}", settings=config["fastp"]
    threads: allocation("fastp", "threads")
    resources:
        mem_mb=allocation("fastp", "mem_mb")
    log: "logs/fastp/{sample}.log"
    benchmark: "benchmarks/fastp/{sample}.tsv"
    conda: "../envs/reads_qc.yaml"
    script: "../scripts/dispatch.py"


if config["fastqc"]["enabled"]:
    rule fastqc_raw:
        input:
            preflight=f"{OUT}/provenance",
            reads=lambda wc: str(getattr(SAMPLES[wc.sample], f"illumina_{wc.mate.lower()}_reads"))
        output:
            directory(f"{OUT}/qc/fastqc/{{sample}}/raw_{{mate}}")
        wildcard_constraints:
            mate="R[12]"
        params:
            stage="fastqc", outdir=lambda wc: f"{OUT}/qc/fastqc/{wc.sample}/raw_{wc.mate}"
        threads: allocation("fastqc", "threads")
        resources:
            mem_mb=allocation("fastqc", "mem_mb")
        log: "logs/fastqc/{sample}.raw.{mate}.log"
        benchmark: "benchmarks/fastqc/{sample}.raw.{mate}.tsv"
        conda: "../envs/reads_qc.yaml"
        script: "../scripts/dispatch.py"

    rule fastqc_clean:
        input:
            reads=f"{OUT}/qc/fastp/{{sample}}/{{mate}}.fastq.gz"
        output:
            directory(f"{OUT}/qc/fastqc/{{sample}}/clean_{{mate}}")
        wildcard_constraints:
            mate="R[12]"
        params:
            stage="fastqc", outdir=lambda wc: f"{OUT}/qc/fastqc/{wc.sample}/clean_{wc.mate}"
        threads: allocation("fastqc", "threads")
        resources:
            mem_mb=allocation("fastqc", "mem_mb")
        log: "logs/fastqc/{sample}.clean.{mate}.log"
        benchmark: "benchmarks/fastqc/{sample}.clean.{mate}.tsv"
        conda: "../envs/reads_qc.yaml"
        script: "../scripts/dispatch.py"


if config["nanoplot"]["enabled"]:
    rule nanoplot:
        input:
            preflight=f"{OUT}/provenance",
            reads=lambda wc: str(SAMPLES[wc.sample].nanopore_reads)
        output:
            directory(f"{OUT}/qc/nanoplot/{{sample}}")
        params:
            stage="nanoplot", outdir=lambda wc: f"{OUT}/qc/nanoplot/{wc.sample}"
        threads: allocation("nanoplot", "threads")
        resources:
            mem_mb=allocation("nanoplot", "mem_mb")
        log: "logs/nanoplot/{sample}.log"
        benchmark: "benchmarks/nanoplot/{sample}.tsv"
        conda: "../envs/reads_qc.yaml"
        script: "../scripts/dispatch.py"


QC_REPORTS = expand(f"{OUT}/qc/fastp/{{sample}}/fastp.json", sample=IDS)
if config["fastqc"]["enabled"]:
    QC_REPORTS += expand(f"{OUT}/qc/fastqc/{{sample}}/{{phase}}_{{mate}}",
                         sample=IDS, phase=["raw", "clean"], mate=["R1", "R2"])
if config["nanoplot"]["enabled"]:
    QC_REPORTS += expand(f"{OUT}/qc/nanoplot/{{sample}}", sample=ONT)

rule multiqc:
    input:
        reports=QC_REPORTS
    output:
        directory(f"{OUT}/qc/multiqc")
    params:
        stage="multiqc", outdir=lambda wc, output: str(output[0])
    threads: allocation("multiqc", "threads")
    resources:
        mem_mb=allocation("multiqc", "mem_mb")
    log: "logs/multiqc.log"
    benchmark: "benchmarks/multiqc.tsv"
    conda: "../envs/reads_qc.yaml"
    script: "../scripts/dispatch.py"

rule reads_qc_all:
    input:
        f"{OUT}/qc/multiqc"
