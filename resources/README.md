# External database setup

These commands are **explicit setup operations**, never workflow rules. They download
large resources: inspect the linked release sizes, available disk and licensing first.
No fixed disk estimate is given because database releases and full/light choices vary.
Use the matching tool environment for every command. Run from the repository root.
Commands below create dedicated named setup environments (optional duplicates of
Snakemake environments); remove them after setup if storage is limited.

## CheckM2 1.1.0

Official [installation](https://github.com/chklovski/CheckM2/tree/1.1.0) and
[compatible reference archive](https://zenodo.org/records/14897628).

```bash
mamba env create -n bactgen-checkm2-setup -f workflow/envs/checkm2.yaml
mkdir -p resources/databases/checkm2
conda run -n bactgen-checkm2-setup checkm2 database --download --path resources/databases/checkm2
conda run -n bactgen-checkm2-setup checkm2 --version
```

Locate the resulting `uniref100.KO.1.dmnd`; set `databases.checkm2.path` to the **file**,
not its directory. CheckM2 1.1.0 validates database compatibility; use the linked
reference release if the current downloader changes. With the resolved file path:

```bash
CHECKM2DB=/absolute/path/to/uniref100.KO.1.dmnd \
conda run -n bactgen-checkm2-setup checkm2 testrun
```

Record the Zenodo record/version and its provided checksum plus a local checksum.
The database is passed via `predict --database_path`; it is never downloaded on demand.

## Bakta 1.10.4

Official [database documentation](https://github.com/oschwengers/bakta/tree/v1.10.4#database)
and [versioned releases](https://doi.org/10.5281/zenodo.4247253).

```bash
mamba env create -n bactgen-bakta-setup -f workflow/envs/bakta.yaml
mkdir -p resources/databases/bakta
conda run -n bactgen-bakta-setup bakta_db list
conda run -n bactgen-bakta-setup bakta_db download --output resources/databases/bakta --type light
conda run -n bactgen-bakta-setup bakta --version
```

The helper selects a compatible schema (Bakta 1.10.4 requires schema 5.x). `--type full` can be chosen for broader
annotation; light/full tradeoffs are documented upstream. Set the actual extracted
directory as `databases.bakta.path`; copy the release/schema/type from its `version.json`
into `.version`. For a repeatable download, use `bakta_db download --minor N --type light`
with the minor release you recorded (the installed Bakta determines the major schema).
Bakta checks database compatibility on invocation. Use the integration test to verify
full runtime dependencies (including external predictors), not just file presence.

## AMRFinderPlus 3.12.8

Official [update instructions](https://github.com/ncbi/amr/wiki/AMRFinderPlus-database)
and [database archives](https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/).

```bash
mamba env create -n bactgen-amrfinder-setup -f workflow/envs/amrfinder.yaml
mkdir -p resources/databases/amrfinder
conda run -n bactgen-amrfinder-setup amrfinder_update --database resources/databases/amrfinder
conda run -n bactgen-amrfinder-setup amrfinder --version
ls resources/databases/amrfinder
```

The updater selects a compatible release and creates a dated version directory.
Point `databases.amrfinder.path` to that **concrete dated directory**, avoiding `latest`
symlinks for published analyses; record its date/version. The directory must contain
indexed BLAST/HMM files, not only the raw reference FASTAs. Confirm metadata with:

```bash
conda run -n bactgen-amrfinder-setup amrfinder --database resources/databases/amrfinder/RELEASE --database_version
conda run -n bactgen-amrfinder-setup amrfinder --database resources/databases/amrfinder/RELEASE --list
```

`RELEASE` is a value to replace with the actual release directory. Use `--list` output
before setting an `amrfinder.organism` value; a genus/species name is not automatically
a supported organism option. Nucleotide mode searches translated DNA without protein
HMM evidence; it is supported and independent of Bakta annotation formatting.

## geNomad 1.8.1

Official [setup](https://portal.nersc.gov/genomad/quickstart.html) and
[database releases](https://portal.nersc.gov/genomad/__data__/) (verify release compatibility
against the installed software; preserve the downloaded release).

```bash
mamba env create -n bactgen-genomad-setup -f workflow/envs/genomad.yaml
mkdir -p resources/databases/genomad
conda run -n bactgen-genomad-setup genomad download-database resources/databases/genomad
conda run -n bactgen-genomad-setup genomad --version
```

Set `databases.genomad.path` to the resulting `genomad_db` directory and record the
version in its metadata/release archive. Store the archive checksum. Run a small
`genomad end-to-end --threads 2 assembly.fasta output/path DB_PATH` to validate that
specific installation if desired. This is an explicit biological analysis, not a
required download step. Consult the pinned software release if a newer database
fails compatibility checks. The workflow never guesses database version from a path.

## ISEScan 1.7.3

Profile models `clusters.faa.hmm` and `clusters.single.faa` are distributed with the
[ISEScan package](https://github.com/xiezhq/ISEScan/tree/v1.7.3). No additional
external database path is required. The isolated environment pins Python 3.9,
NumPy 1.23.5, SciPy 1.10.1, pandas 1.5.3 and fastcluster 1.2.6 to avoid NumPy 2 ABI
and removed-API problems in this older program. Environment resolution remains
platform-dependent and has not been verified here.

```bash
mamba env create -n bactgen-isescan-setup -f workflow/envs/isescan.yaml
conda run -n bactgen-isescan-setup isescan.py --help
conda run -n bactgen-isescan-setup hmmsearch -h
conda run -n bactgen-isescan-setup phmmer -h
```

ISEScan also uses FragGeneScan and BLAST. A real integration run, including both
completed HMMER search outputs, is required to establish working installation.
Do not interpret a zero process return code alone as evidence of a valid no-hit result.

## Freeze and update

1. Download into a new release-specific directory. Keep any previous release intact.
2. Record upstream release identifiers, download URLs and archive SHA-256 checksums.
3. Write a checksum inventory of the installed database if practical; hashing large
   databases is a separate potentially lengthy operation, never done automatically.
4. Update `config.yaml` paths and version strings; run structural and database validation.
5. Resolve/export each environment for your platform, for example
   `conda list --explicit -n bactgen-isescan-setup > isescan-linux-64.explicit.txt`.
6. Run the real-data integration test and review changed calls before using new releases.

Preflight verifies version declarations and database existence/non-emptiness; it does
not verify checksums, archive identity or all internal indices. Biological tools do
additional native compatibility checks. Keep databases read-only while jobs run.
