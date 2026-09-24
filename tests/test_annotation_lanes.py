"""Annotation-lane results: InterPro domain sharing, BioCyc/MetaCyc
co-membership (honest weak signal), Rhea design reactions, OLS4 ChEBI
terms, Crossref bibliography verification."""
import csv
import json
from pathlib import Path

RES = Path(__file__).resolve().parent.parent / "results"


def test_interpro_domain_sharing_strong():
    d = json.loads((RES / "interpro_domain_sharing.json").read_text())
    assert d["n_mapped"] >= 20
    assert d["pairs_sharing_domain"] >= 10
    assert d["fisher_p"] < 1e-10


def test_metacyc_co_membership_honest_weak():
    d = json.loads((RES / "metacyc_pathway_summary.json").read_text())
    assert d["n_pairs"] == 14
    assert d["pairs_co_member"] <= 2  # sparse BioCyc yeast coverage: no strong claim


def test_rhea_design_reactions():
    rows = list(csv.DictReader(open(RES / "rhea_design_reactions.csv")))
    assert len(rows) >= 6
    ecs = {r["ec"] for r in rows}
    assert any("1.3.5.1" in e for e in ecs) and any("3.1.2.1" in e for e in ecs)


def test_ols_chebi_terms():
    rows = list(csv.DictReader(open(RES / "ols_chebi_design_metabolites.csv")))
    assert len(rows) == 5
    by_q = {r["query"]: r["chebi_id"] for r in rows}
    assert by_q["succinate"] == "CHEBI:26806"
    assert all(v.startswith("CHEBI:") for v in by_q.values())


def test_crossref_bibliography():
    d = json.loads((RES / "crossref_bib_summary.json").read_text())
    assert d["n_references"] == 14
    assert d["n_crossref_matched"] >= 12


def test_orthodb_same_group_replicates():
    d = json.loads((RES / "orthodb_same_group.json").read_text())
    assert d["n_mapped"] >= 20
    assert d["pairs_same_group"] >= 10
    assert d["fisher_p"] < 1e-8
