# Optional real-data integration

Run on Linux x86-64 after installing the controller environment and compatible
biological databases. This is not part of lightweight CI.

Unicycler distributes an explicitly small example consisting of synthetic
short and long reads from Shigella sonnei plasmids A, B and E. These are synthetic
example reads, not experimental Nanopore data. The assembler accepts them as
long-read FASTQ; this checks hybrid wiring, not platform-specific Nanopore accuracy.
See the [original data description](https://github.com/rrwick/Unicycler/blob/v0.5.1/sample_data/README.md).
For a Nanopore-specific validation, supply a real matched Nanopore/Illumina isolate.

Explicitly download the small version-pinned example (about 7.7 MB total):

```bash
mkdir -p data/raw/integration
curl -fL https://raw.githubusercontent.com/rrwick/Unicycler/v0.5.1/sample_data/short_reads_1.fastq.gz -o data/raw/integration/R1.fastq.gz
curl -fL https://raw.githubusercontent.com/rrwick/Unicycler/v0.5.1/sample_data/short_reads_2.fastq.gz -o data/raw/integration/R2.fastq.gz
curl -fL https://raw.githubusercontent.com/rrwick/Unicycler/v0.5.1/sample_data/long_reads_low_depth.fastq.gz -o data/raw/integration/long.fastq.gz
sha256sum data/raw/integration/*.gz > data/raw/integration/reads.sha256
```

Create `manifest.integration.csv` in the repository root:

```csv
sample,nanopore_reads,illumina_r1_reads,illumina_r2_reads
integration,data/raw/integration/long.fastq.gz,data/raw/integration/R1.fastq.gz,data/raw/integration/R2.fastq.gz
```

Copy `config.yaml` to `config.integration.yaml`, set `manifest: manifest.integration.csv`
and `results_dir: results/integration`, and fill all database paths/versions. Leave
`checkm2.exclude_failed: false`: these plasmids are not a complete bacterial genome.
Completeness is expected to be low and is not an accuracy benchmark on this input.
Leave all modules enabled to exercise the entire graph. Also repeat with a real
Illumina-only sample by leaving its Nanopore column empty and using a new results path.

```bash
python workflow/scripts/validate_manifest.py manifest.integration.csv
python workflow/scripts/validate_config.py config.integration.yaml --check-databases
snakemake --configfile config.integration.yaml --cores 4 --dry-run
BACTGEN_INTEGRATION_CONFIG="$PWD/config.integration.yaml" BACTGEN_CORES=4 pytest -v -m integration
```

Inspect assembly.fasta, every native output and all summary statuses. Compare a
known complete bacterial isolate against an independently run reference workflow
before drawing biological conclusions. This project has not run that validation
in the authoring environment. Real analysis may take substantial time despite small
input reads, because several tools load large databases/models.
