"""Reproduce the isozyme over-rescue finding on the in-repo data."""
import pandas as pd

from yeasttwin.labels import load_labels
from yeasttwin.overrescue import overrescue_audit


def test_overrescue_finding_reproduces():
    ko = pd.Series(pd.read_json("results/fba_single_ko.json")["ratio"])
    out = overrescue_audit(ko, load_labels())
    assert out["n_essential"] == 159
    assert out["missed_isozyme_backed_frac"] > 0.15
    assert out["caught_isozyme_backed_frac"] < 0.05
    assert out["fisher_odds_ratio"] > 5
    assert out["fisher_p"] < 0.01
    assert len(out["overrescued_genes"]) >= 15
