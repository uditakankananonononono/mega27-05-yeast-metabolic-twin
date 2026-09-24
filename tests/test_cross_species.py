"""Hermetic tests for the cross-species over-rescue audit."""
import json
from pathlib import Path

import pandas as pd

from yeasttwin.cross_species import load_keio_essentials

ROOT = Path(__file__).resolve().parents[1]


def test_keio_table_parses():
    b = load_keio_essentials()
    assert 295 <= len(b) <= 305  # ~303 essential candidates
    assert all(x.startswith("b") and len(x) == 5 for x in b)
    assert "b0048" in b  # folA, canonical essential


def test_ecoli_audit_result_is_significant_and_consistent():
    out = json.loads((ROOT / "results" / "overrescue_audit_ecoli.json")
                     .read_text())
    assert out["n_essential"] == out["n_missed"] + out["n_caught"]
    assert out["fisher_p"] < 0.001
    assert out["missed_isozyme_backed_frac"] > \
        5 * out["caught_isozyme_backed_frac"]
    for g in out["overrescued_genes"]:
        assert g.startswith("b")
