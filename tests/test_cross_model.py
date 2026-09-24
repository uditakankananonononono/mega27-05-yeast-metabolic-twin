"""Hermetic test: cross-reconstruction replication (iMM904)."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_imm904_replication():
    out = json.loads((ROOT / "results" / "overrescue_audit_imm904.json")
                     .read_text())
    assert out["fisher_p"] < 0.01
    assert out["caught_isozyme_backed_frac"] == 0.0
    gem = set(json.loads((ROOT / "results" / "overrescue_audit.json")
                         .read_text())["overrescued_genes"])
    overlap = gem & set(out["overrescued_genes"])
    assert len(overlap) >= 10  # same genes recur across reconstructions
