from pathlib import Path
from types import SimpleNamespace

import pytest
from conftest import ROOT
from run_stage import execute


def job(tmp_path, config, stage):
    return SimpleNamespace(
        scriptdir=str(ROOT / "workflow/scripts"),
        config=config,
        params=SimpleNamespace(stage=stage, outdir=str(tmp_path / "output")),
        wildcards=SimpleNamespace(sample="S1"),
        threads=1,
        log=[str(tmp_path / "stage.log")],
        input=SimpleNamespace(assembly="test.fasta"),
    )


def test_failed_process_never_promoted(tmp_path, config, monkeypatch):
    import subprocess

    def fail(*args, **kwargs):
        raise subprocess.CalledProcessError(1, args[0])

    monkeypatch.setattr(subprocess, "run", fail)
    smk = job(tmp_path, config, "amrfinder")
    with pytest.raises(subprocess.CalledProcessError):
        execute(smk)
    assert not Path(smk.params.outdir).exists()
    assert "CalledProcessError" in Path(smk.log[0]).read_text()


def test_success_without_required_output_rejected(tmp_path, config, monkeypatch):
    import subprocess

    monkeypatch.setattr(subprocess, "run", lambda *a, **k: None)
    smk = job(tmp_path, config, "amrfinder")
    with pytest.raises(FileNotFoundError):
        execute(smk)
    assert not Path(smk.params.outdir).exists()
