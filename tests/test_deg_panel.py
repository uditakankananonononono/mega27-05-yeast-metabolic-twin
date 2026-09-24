"""DEG panel: 3 new species + 4 cross-study replications (committed results)."""
import json
from pathlib import Path

import pandas as pd

RESULTS = Path(__file__).resolve().parents[1] / "results"


def load(suffix):
    return json.loads((RESULTS / f"overrescue_audit_{suffix}.json").read_text())


def test_new_species_significant():
    """S. aureus (Ji 2001) and S. Typhimurium (Knuth 2004) replicate."""
    sau = load("saureus_deg1002")
    assert sau["fisher_p"] < 0.05 and len(sau["overrescued_genes"]) == 10
    assert all(g.startswith("SA_RS") for g in sau["overrescued_genes"])
    stm = load("styphimurium_deg1011")
    assert stm["fisher_p"] < 1e-3 and len(stm["overrescued_genes"]) == 14
    assert all(g.startswith("STM") for g in stm["overrescued_genes"])


def test_cross_study_replications():
    """Independent essentiality studies replicate in E. coli, B. subtilis, Mtb."""
    ec = load("ecoli_gerdes_deg1018")       # Gerdes 2003, independent of Keio
    assert 5 < ec["fisher_odds_ratio"] and ec["fisher_p"] < 1e-5
    assert len(ec["overrescued_genes"]) == 45
    bs = load("bsubtilis_kobayashi_deg1001")  # Kobayashi 2003, independent of Koo
    assert bs["fisher_p"] < 5e-3 and len(bs["overrescued_genes"]) == 8
    mtb1 = load("mtb_griffin_deg1025")      # Griffin 2011 Tn-seq, independent
    assert mtb1["fisher_p"] < 1e-6 and len(mtb1["overrescued_genes"]) == 48
    mtb2 = load("mtb_zhang_deg1027")        # Zhang 2012, independent
    assert mtb2["fisher_p"] < 1e-5 and len(mtb2["overrescued_genes"]) == 44


def test_hpylori_honest_underpowered():
    """H. pylori iIT341 is too small to test - recorded honestly, not spun."""
    hp = load("hpylori_deg1008")
    assert hp["n_in_model"] <= 40
    assert hp["fisher_p"] > 0.05
    assert len(hp["overrescued_genes"]) <= 2


def test_summary_table():
    df = pd.read_csv(RESULTS / "deg_panel_summary.csv")
    assert len(df) == 7
    sig = df[df.fisher_p < 0.05]
    assert len(sig) == 6  # all except H. pylori
