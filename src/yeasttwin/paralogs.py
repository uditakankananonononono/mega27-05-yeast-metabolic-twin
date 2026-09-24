"""Mechanism analysis for isozyme over-rescue: sequence identity between
each over-rescued essential gene and its model 'backup' isozyme partners.

High identity => annotation-driven isozyme (paralog transfer) - the likely
false-backup mechanism. Low identity => a genuinely different enzyme the
model trusts to catalyse the same reaction.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import pandas as pd
from Bio import Align

from .evaluate import load_sequences
from .model import load_model

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results"


def isozyme_partners(model, gene_id: str) -> set[str]:
    """Genes appearing in ALTERNATIVE GPR clauses of the gene's reactions."""
    partners: set[str] = set()
    for rxn in model.genes.get_by_id(gene_id).reactions:
        rule = rxn.gene_reaction_rule.strip()
        if not rule or " or " not in rule:
            continue
        clauses = re.split(r"\bor\b", rule)
        own = [c for c in clauses if gene_id in c]
        others = [c for c in clauses if gene_id not in c]
        for c in others:
            for g in re.findall(r"[A-Za-z0-9_\-]+", c.replace("and", " ")):
                if g in model.genes and g != gene_id:
                    partners.add(g)
    return partners


def identity(a: str, b: str) -> float | None:
    if not a or not b:
        return None
    aln = Align.PairwiseAligner()
    aln.mode = "global"
    aln.match_score = 1.0
    aln.mismatch_score = 0.0
    aln.open_gap_score = -0.5
    aln.extend_gap_score = -0.1
    score = aln.score(a, b)
    return float(score / max(len(a), len(b)))


def main():
    model = load_model()
    genes = pd.read_csv(RESULTS / "overrescued_genes.csv")["gene"]
    seqs = dict(zip([g.id for g in model.genes],
                    load_sequences([g.id for g in model.genes])))
    rows = []
    for gid in genes:
        partners = sorted(isozyme_partners(model, gid))
        idents = [identity(seqs.get(gid, ""), seqs.get(p, "")) for p in partners]
        idents = [x for x in idents if x is not None]
        rows.append({
            "gene": gid,
            "n_partners": len(partners),
            "partners": ";".join(partners),
            "max_identity": max(idents) if idents else None,
            "mean_identity": sum(idents) / len(idents) if idents else None,
        })
    df = pd.DataFrame(rows)
    df.to_csv(RESULTS / "overrescued_paralog_identity.csv", index=False)
    print(df.to_string(index=False))
    hi = df[df.max_identity > 0.3]
    print("\n%d/%d over-rescued genes have a close-paralog 'backup' (identity>0.30)"
          % (len(hi), len(df)))


if __name__ == "__main__":
    main()
