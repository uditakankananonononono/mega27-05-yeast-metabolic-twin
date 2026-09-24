"""Hermetic test: B. subtilis abundance honest negative."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_bsub_abundance_negative():
    out = json.loads((ROOT / "results" / "overrescue_abundance_bsubtilis.json")
                     .read_text())
    assert out["n_pairs"] >= 60
    assert 0.4 < out["frac_backup_lower"] < 0.7  # NOT replicated
    assert out["wilcoxon_p"] > 0.05  # honestly non-significant
