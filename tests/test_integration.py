"""Opt-in real workflow execution; never uses unit-test reads as biological inputs."""

import os
import subprocess
import sys
from pathlib import Path

import pytest
from conftest import ROOT


@pytest.mark.integration
@pytest.mark.skipif(
    not os.environ.get("BACTGEN_INTEGRATION_CONFIG"),
    reason="Set BACTGEN_INTEGRATION_CONFIG to a prepared real-data configuration",
)
def test_real_pipeline():
    config = Path(os.environ["BACTGEN_INTEGRATION_CONFIG"]).resolve()
    subprocess.run(
        [
            sys.executable,
            "-m",
            "snakemake",
            "--configfile",
            str(config),
            "--use-conda",
            "--cores",
            os.environ.get("BACTGEN_CORES", "4"),
            "--rerun-incomplete",
        ],
        cwd=ROOT,
        check=True,
    )
    second = subprocess.run(
        [
            sys.executable,
            "-m",
            "snakemake",
            "--configfile",
            str(config),
            "--cores",
            "1",
            "--dry-run",
        ],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    assert "Nothing to be done" in second.stdout
