"""Hermetic test: SGD inviable triangulation result."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_triangulation():
    out = json.loads((ROOT / "results" / "sgd_inviable_triangulation.json")
                     .read_text())
    assert out["n_overrescued"] == 19
    assert out["n_inviable_in_sgd"] == 18
    assert out["missing"] == ["YMR108W"]  # ILV2 curation conflict
