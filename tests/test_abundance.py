"""Hermetic test for the PaxDb abundance mechanism analysis."""
import json
from pathlib import Path

import pandas as pd

from yeasttwin.abundance import load_paxdb

ROOT = Path(__file__).resolve().parents[1]


def test_paxdb_loads():
    ppm = load_paxdb()
    assert len(ppm) > 5000
    assert ppm["YGL008C"] > ppm["YER005W"]  # PMA1 >> backup


def test_abundance_result_supports_low_expression_backups():
    out = json.loads((ROOT / "results" / "overrescue_abundance.json").read_text())
    assert out["n_pairs"] >= 25
    assert out["frac_backup_lower"] > 0.5
    assert out["wilcoxon_p"] < 0.05
