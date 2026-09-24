"""Per-gene mechanistic features derived from the genome-scale model.

Every feature is computed from model structure and constraint-based
simulation only; experimental essentiality labels are never consulted, so
features cannot leak labels.
"""
from __future__ import annotations

import re
from collections import defaultdict

import cobra
import numpy as np
import pandas as pd
from cobra.flux_analysis import pfba, single_gene_deletion

# Highly connected cofactors that create spurious shortcuts in metabolite graphs.
CURRENCY_NAMES = {
    "H2O", "H+", "ATP", "ADP", "AMP", "NAD", "NADH", "NADP(+)", "NADPH",
    "phosphate", "diphosphate", "carbon dioxide", "oxygen", "coenzyme A",
    "ammonium", "FAD", "FADH2", "NADP", "NAD(+)", "NADP(H)",
}


def is_currency(met: cobra.Metabolite) -> bool:
    return met.name in CURRENCY_NAMES


def gpr_stats(model: cobra.Model) -> pd.DataFrame:
    """Isozyme redundancy and complex membership for each gene.

    For each reaction we parse the GPR into OR-clauses of AND-complexes.
    ``min_isozymes`` = fewest alternative clauses among the gene's reactions
    (1 means the gene sits in a reaction with no backup enzyme).
    """
    rows = defaultdict(lambda: {"n_rxns": 0, "min_isozymes": 99, "max_complex": 0,
                                "n_sole_rxns": 0})
    for rxn in model.reactions:
        rule = rxn.gene_reaction_rule.strip()
        if not rule:
            continue
        clauses = [c for c in re.split(r"\bor\b", rule)]
        n_alt = len(clauses)
        for clause in clauses:
            genes = re.findall(r"[A-Za-z0-9_\-]+", clause.replace("and", " "))
            genes = [g for g in genes if g in model.genes]
            for g in genes:
                r = rows[g]
                r["n_rxns"] += 1
                r["min_isozymes"] = min(r["min_isozymes"], n_alt)
                r["max_complex"] = max(r["max_complex"], len(genes))
                if n_alt == 1:
                    r["n_sole_rxns"] += 1
    df = pd.DataFrame.from_dict(rows, orient="index")
    df.index.name = "gene"
    return df


def flux_features(model: cobra.Model) -> pd.DataFrame:
    """pFBA wild-type flux carried by each gene's reactions."""
    sol = pfba(model)
    out = {}
    for g in model.genes:
        fl = [abs(sol.fluxes[r.id]) for r in g.reactions]
        out[g.id] = {
            "max_abs_flux": max(fl) if fl else 0.0,
            "sum_abs_flux": float(np.sum(fl)) if fl else 0.0,
            "active_rxn_frac": float(np.mean([f > 1e-9 for f in fl])) if fl else 0.0,
        }
    df = pd.DataFrame.from_dict(out, orient="index")
    df.index.name = "gene"
    return df


def knockout_ratios(model: cobra.Model) -> pd.Series:
    """Growth ratio (KO / WT) for every single-gene deletion."""
    wt = model.optimize().objective_value
    ko = single_gene_deletion(model, processes=1)
    ratio = {}
    for _, row in ko.iterrows():
        (gid,) = tuple(row["ids"])
        v = row["growth"]
        v = 0.0 if v != v else float(v)
        ratio[gid] = v / wt if wt > 0 else 0.0
    for g in model.genes:
        ratio.setdefault(g.id, 1.0)
    return pd.Series(ratio, name="ko_ratio")


def topology_features(model: cobra.Model) -> pd.DataFrame:
    """Degree-style features on the gene graph (shared non-currency metabolites)."""
    adj = gene_adjacency(model)
    deg = {g: len(n) for g, n in adj.items()}
    comp = {g.id: len({m.compartment for r in g.reactions for m in r.metabolites})
            for g in model.genes}
    subsys = {g.id: len({(r.subsystem or "") for r in g.reactions}) for g in model.genes}
    df = pd.DataFrame({"gene_degree": pd.Series(deg), "n_compartments": pd.Series(comp),
                       "n_subsystems": pd.Series(subsys)})
    df.index.name = "gene"
    return df.fillna(0)


def gene_adjacency(model: cobra.Model) -> dict[str, set[str]]:
    """Genes are adjacent if their reactions share a non-currency metabolite."""
    met_genes: dict[str, set[str]] = defaultdict(set)
    for rxn in model.reactions:
        gs = {g.id for g in rxn.genes}
        if not gs:
            continue
        for met in rxn.metabolites:
            if not is_currency(met):
                met_genes[met.id] |= gs
    adj: dict[str, set[str]] = {g.id: set() for g in model.genes}
    for gs in met_genes.values():
        if len(gs) > 60:  # hub metabolite; skip to keep graph informative
            continue
        for g in gs:
            adj[g] |= gs - {g}
    return adj
