"""Fix-wave lanes: statsmodels FDR, QuickGO, Ensembl, Complex Portal,
IntAct, PANTHER, RDKit, PRIDE provenance, and the accession-level ledger.
All hermetic - they assert on committed result files only."""
import json
from pathlib import Path

import pandas as pd

RESULTS = Path(__file__).resolve().parents[1] / "results"


def test_statsmodels_fdr_panel():
    df = pd.read_csv(RESULTS / "statsmodels_fdr.csv")
    assert len(df) == 12
    s = json.loads((RESULTS / "statsmodels_fdr_summary.json").read_text())
    assert s["n_significant_fdr05"] >= 10
    assert s["yeast_discovery_padj"] < 1e-3


def test_quickgo_annotations():
    df = pd.read_csv(RESULTS / "quickgo_overrescued.csv")
    assert df.gene.nunique() == 16  # 3 merged/dubious ORFs have no UniProt accession - honest coverage
    assert (df.n_mf_terms > 0).all()
    s = json.loads((RESULTS / "quickgo_summary.json").read_text())
    assert s["frac_with_activity_term"] == 1.0


def test_ensembl_xref_honest_misses():
    df = pd.read_csv(RESULTS / "ensembl_xref_overrescued.csv")
    assert len(df) == 19
    s = json.loads((RESULTS / "ensembl_xref_summary.json").read_text())
    assert s["n_resolved_protein_coding"] == 14
    assert s["n_queried"] - s["n_resolved_protein_coding"] == 5  # recorded


def test_complexportal_intact_pairs():
    cp = pd.read_csv(RESULTS / "complexportal_pairs.csv")
    assert len(cp) == 14
    ia = pd.read_csv(RESULTS / "intact_pairs.csv")
    assert len(ia) == 14
    assert (ia.intact_physical_interactions >= 0).all()
    assert ia.intact_physical_interactions.sum() > 0


def test_panther_rdkit_pride():
    pf = pd.read_csv(RESULTS / "panther_families.csv")
    assert len(pf) == 19 and (pf.family_id.str.len() > 0).all()
    rk = pd.read_csv(RESULTS / "rdkit_design_metabolites.csv")
    assert len(rk) == 5 and rk.tpsa_agree.all()
    assert set(rk["query"][rk.rdkit_formal_charge <= -2]) == {
        "succinate", "fumarate", "alpha-ketoglutarate"}
    pr = pd.read_csv(RESULTS / "pride_yeast_projects.csv")
    s = json.loads((RESULTS / "pride_provenance_summary.json").read_text())
    assert s["n_distinct_pxd"] >= 5
    assert s["n_tokens_with_pride_hits"] >= 10


def test_accession_ledger_gate():
    df = pd.read_csv(RESULTS / "accession_ledger.csv")
    assert not df.duplicated(subset=["resource", "accession"]).any()
    s = json.loads((RESULTS / "accession_ledger_summary.json").read_text())
    assert s["total_distinct_accession_level_records"] >= 120
    assert s["total_distinct_accession_level_records"] == len(df)
    assert s["gate_120"] is True
