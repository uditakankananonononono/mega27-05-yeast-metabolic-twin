"""Hermetic test: E. coli backup abundance replication."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_ecoli_abundance():
    out = json.loads((ROOT / "results" / "overrescue_abundance_ecoli.json")
                     .read_text())
    assert out["n_pairs"] >= 15
    assert out["frac_backup_lower"] > 0.7
    assert out["wilcoxon_p"] < 0.01
    assert out["median_log2_ratio"] > 2.0
