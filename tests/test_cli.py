"""Hermetic test: the overrescue-audit CLI runs end to end."""
import json
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_cli_end_to_end(tmp_path):
    ess = tmp_path / "ess.txt"
    # small essential subset for a fast screen: known over-rescued + caught
    ess.write_text("\n".join(["b0048", "b0126", "b0171", "b0475",
                              "b0733", "b2411", "b2563", "b3640"]) + "\n")
    out = tmp_path / "report.json"
    subprocess.run(
        ["python3", "-m", "yeasttwin.cli", "--model",
         str(ROOT / "data" / "raw" / "iML1515.json"),
         "--essentials", str(ess), "--out", str(out)],
        cwd=ROOT, env={"PYTHONPATH": "src", "PATH": "/usr/bin:/bin"},
        check=True, capture_output=True, timeout=300)
    rep = json.loads(out.read_text())
    assert rep["n_essential"] == 8
    assert rep["n_missed"] + rep["n_caught"] == 8
    assert "fisher_odds_ratio" in rep
