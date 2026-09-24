"""Hermetic test: B. subtilis third-species replication."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_bsub_audit():
    out = json.loads((ROOT / "results" / "overrescue_audit_bsubtilis.json")
                     .read_text())
    assert out["n_essential"] == out["n_missed"] + out["n_caught"]
    assert out["fisher_p"] < 1e-5
    assert out["fisher_odds_ratio"] > 10
    assert out["missed_isozyme_backed_frac"] >= 0.5
