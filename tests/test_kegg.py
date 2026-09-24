"""Hermetic test: KEGG provenance analysis."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_kegg_provenance():
    out = json.loads((ROOT / "results" / "kegg_same_ko_backups.json").read_text())
    assert out["n_genes"] == 19
    assert out["n_same_ko_pairs"] == 14
    joint = json.loads((ROOT / "results" / "kegg_provenance_joint.json").read_text())
    assert joint["sameko_vs_closeparalog_fisher_p"] < 0.001
