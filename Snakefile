from pathlib import Path
import sys

configfile: "config.yaml"
sys.path.insert(0, str(Path(workflow.basedir) / "workflow" / "scripts"))
from validate_config import validate_config
from validate_manifest import validate_manifest

validate_config(config)
SAMPLES = validate_manifest(config["manifest"])
IDS = list(SAMPLES)
ONT = [name for name, record in SAMPLES.items() if record.nanopore_reads]
OUT = config["results_dir"].rstrip("/")
SUMMARY_NAMES = ["assembly_summary", "annotation_summary", "resistance_summary",
                 "mobile_elements_summary", "samples_summary", "resistance_genes",
                 "mobile_elements", "arg_mobile_associations", "checkm2_normalized"]
FINAL = [f"{OUT}/summary/{name}.tsv" for name in SUMMARY_NAMES]

wildcard_constraints:
    sample="[A-Za-z0-9][A-Za-z0-9_.-]*"

rule all:
    input:
        FINAL,
        f"{OUT}/qc/multiqc"

include: "workflow/rules/common.smk"
include: "workflow/rules/preflight.smk"
include: "workflow/rules/reads_qc.smk"
include: "workflow/rules/assembly.smk"
include: "workflow/rules/assembly_qc.smk"
include: "workflow/rules/checkm2.smk"
include: "workflow/rules/bakta.smk"
include: "workflow/rules/amrfinder.smk"
include: "workflow/rules/genomad.smk"
include: "workflow/rules/isescan.smk"
include: "workflow/rules/summary.smk"
