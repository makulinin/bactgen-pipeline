from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "tests" / "data"


@pytest.fixture
def config():
    return yaml.safe_load((ROOT / "config.yaml").read_text())


@pytest.fixture
def result_tree(tmp_path, config):
    import json
    import shutil

    config["results_dir"] = str(tmp_path / "results")
    root = Path(config["results_dir"])
    paths = {
        "quast/S1/report.tsv": "quast.tsv",
        "checkm2/quality_report.tsv": "checkm2.tsv",
        "bakta/S1/annotation.gff3": "annotation.gff3",
        "amrfinder/S1/amrfinder.tsv": "amrfinder.tsv",
        "genomad/S1/assembly_summary/assembly_plasmid_summary.tsv": "plasmid.tsv",
        "genomad/S1/assembly_summary/assembly_virus_summary.tsv": "virus.tsv",
        "isescan/S1/prediction/assembly.fasta.csv": "isescan.csv",
    }
    for target, source in paths.items():
        file = root / target
        file.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(DATA / source, file)
    fasta = root / "assembly/S1/assembly.fasta"
    fasta.parent.mkdir(parents=True)
    fasta.write_text(">c1\n" + "ACGT" * 300 + "\n>c2\n" + "ACGT" * 200 + "\n")
    for module in ("bakta", "amrfinder", "genomad", "isescan"):
        state = {"status": "success"}
        if module == "isescan":
            state["native_table"] = "prediction/assembly.fasta.csv"
        (root / module / "S1/status.json").write_text(json.dumps(state))
    return root, config
