"""Hermetic tests: M. tuberculosis audit + E. coli K-12 panel."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_mtb_audit():
    out = json.loads((ROOT / "results" / "overrescue_audit_mtb.json").read_text())
    assert out["n_essential"] == out["n_missed"] + out["n_caught"]
    assert out["fisher_p"] < 1e-4
    assert out["fisher_odds_ratio"] > 5
    assert len(out["overrescued_genes"]) >= 30


def test_k12_panel():
    s = json.loads((ROOT / "results" / "ecoli_strain_panel_summary.json")
                   .read_text())
    assert s["models_audited"] == 4
    assert s["models_significant"] == 3
    assert s["recurrent_genes"]["b0733"] == 4  # recurs in every generation
