"""Snakemake script entry point, compatible with injected Snakemake preambles."""

import sys

sys.path.insert(0, snakemake.scriptdir)  # noqa: F821
from run_stage import execute  # noqa: E402

execute(snakemake)  # noqa: F821
