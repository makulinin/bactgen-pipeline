import subprocess
import sys

import pytest
import yaml
from conftest import ROOT


def snakemake(*args):
    return subprocess.run(
        [sys.executable, "-m", "snakemake", "--cores", "2", *args],
        cwd=ROOT,
        capture_output=True,
        text=True,
        timeout=90,
    )


def test_full_dag():
    result = snakemake("--configfile", "tests/config.yaml", "--dry-run")
    assert result.returncode == 0, result.stdout + result.stderr
    for stage in (
        "assemble",
        "fastp",
        "quast",
        "checkm2",
        "bakta",
        "amrfinder",
        "genomad",
        "isescan",
        "nanoplot",
        "multiqc",
        "summarize",
    ):
        assert stage in result.stdout
    graph = snakemake("--configfile", "tests/config.yaml", "--dag")
    assert graph.returncode == 0, graph.stderr
    assert "digraph snakemake_dag" in graph.stdout


@pytest.mark.parametrize(
    "target", ["reads_qc_all", "assembly_all", "assembly_qc_all", "annotation_all", "summary_all"]
)
def test_named_targets(target):
    result = snakemake(target, "--configfile", "tests/config.yaml", "--dry-run")
    assert result.returncode == 0, result.stdout + result.stderr


@pytest.mark.parametrize("enabled,exclude", [(False, False), (True, True)])
def test_optional_and_quality_gate(tmp_path, enabled, exclude):
    cfg = yaml.safe_load((ROOT / "tests/config.yaml").read_text())
    cfg["modules"] = dict.fromkeys(cfg["modules"], enabled)
    cfg["checkm2"]["exclude_failed"] = exclude
    cfg["fastqc"]["enabled"] = enabled
    cfg["nanoplot"]["enabled"] = enabled
    cfg["results_dir"] = "custom-results"
    file = tmp_path / "config.yaml"
    file.write_text(yaml.safe_dump(cfg))
    result = snakemake("--configfile", str(file), "--dry-run")
    assert result.returncode == 0, result.stdout + result.stderr
    if not enabled:
        assert "rule bakta:" not in result.stdout
        assert "rule checkm2:" not in result.stdout
    else:
        assert "custom-results/checkm2" in result.stdout


def test_real_preflight_dispatch_and_reuse(tmp_path):
    output = tmp_path / "provenance"
    log = tmp_path / "preflight.log"
    config = yaml.safe_load((ROOT / "tests/config.yaml").read_text())
    config["modules"] = dict.fromkeys(config["modules"], False)
    cfg_file = tmp_path / "config.yaml"
    cfg_file.write_text(yaml.safe_dump(config))
    snakefile = tmp_path / "Snakefile"
    snakefile.write_text(f"""
rule preflight:
    input: "tests/data/manifest.csv"
    output: directory({str(output)!r})
    params: stage="preflight", outdir={str(output)!r}
    threads: 1
    resources: mem_mb=1024
    log: {str(log)!r}
    script: {str(ROOT / "workflow/scripts/dispatch.py")!r}
""")
    # This tests the actual script preamble without installing biological environments.
    result = snakemake("--snakefile", str(snakefile), "--configfile", str(cfg_file))
    assert result.returncode == 0, result.stdout + result.stderr
    assert (output / "resolved_manifest.json").is_file()
    assert (output / "status.json").is_file()
    repeat = snakemake("--snakefile", str(snakefile), "--configfile", str(cfg_file), "--dry-run")
    assert repeat.returncode == 0, repeat.stdout + repeat.stderr
    assert "Nothing to be done" in repeat.stdout


@pytest.mark.parametrize("module", ["checkm2", "bakta", "amrfinder", "genomad", "isescan"])
def test_independent_module_disable(tmp_path, module):
    config = yaml.safe_load((ROOT / "tests/config.yaml").read_text())
    config["modules"][module] = False
    file = tmp_path / "config.yaml"
    file.write_text(yaml.safe_dump(config))
    result = snakemake("--configfile", str(file), "--dry-run")
    assert result.returncode == 0, result.stdout + result.stderr
    assert f"rule {module}:" not in result.stdout
