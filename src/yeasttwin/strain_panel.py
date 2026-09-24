"""E. coli K-12 reconstruction panel: over-rescue across model generations.

Every BiGG E. coli strain reconstruction whose genes use b-numbers is
audited against the same Keio essential-gene list. Question: is the
over-rescue signal stable across independent reconstructions of one
species, and do the same genes recur?
"""
from __future__ import annotations

import json
import re
import urllib.request
from pathlib import Path

import cobra
import pandas as pd
from cobra.flux_analysis import single_gene_deletion

from .cross_species import load_keio_essentials
from .overrescue import overrescue_audit

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data" / "raw"
RESULTS = ROOT / "results"
MODELS_DIR = DATA / "ecoli_strain_models"

PANEL = ["e_coli_core", "iAF1260", "iAF1260b", "iJO1366",
         "iECDH10B_1368", "iECDH1ME8569_1439", "iEcolC_1368"]


def fetch(mid: str) -> Path:
    MODELS_DIR.mkdir(exist_ok=True)
    p = MODELS_DIR / f"{mid}.json"
    if not p.exists():
        url = f"http://bigg.ucsd.edu/static/models/{mid}.json"
        with urllib.request.urlopen(url, timeout=60) as r:
            p.write_bytes(r.read())
    return p


def bnum_frac(model) -> float:
    ids = [g.id for g in model.genes]
    return sum(bool(re.fullmatch(r"b\d{4}", g)) for g in ids) / max(len(ids), 1)


def audit_one(mid: str, keio: list[str]) -> dict | None:
    try:
        model = cobra.io.load_json_model(str(fetch(mid)))
    except Exception as e:
        return {"model": mid, "error": str(e)[:80]}
    if bnum_frac(model) < 0.5:
        return {"model": mid, "skipped": "gene ids not b-numbers"}
    wt = model.slim_optimize()
    if not wt or wt < 1e-6:
        return {"model": mid, "skipped": f"no growth (WT={wt})"}
    ess = [b for b in keio if b in model.genes]
    ko = single_gene_deletion(model, gene_list=ess, processes=1)
    ko = ko.dropna(subset=["growth"])
    ratios = pd.Series({next(iter(r["ids"])): float(r["growth"]) / wt
                        for _, r in ko.iterrows()})
    out = overrescue_audit(ratios, pd.Series(True, index=ess), model=model)
    out["model"] = mid
    out["wt"] = float(wt)
    return out


def main():
    keio = load_keio_essentials()
    rows = []
    for i, mid in enumerate(PANEL):
        r = audit_one(mid, keio)
        rows.append(r)
        tag = (f"OR {r['fisher_odds_ratio']:.1f} p {r['fisher_p']:.1e}"
               if r and "fisher_odds_ratio" in r
               else str(r)[:60])
        print(f"[{i+1}/{len(PANEL)}] {mid}: {tag}", flush=True)
        pd.DataFrame(rows).to_csv(RESULTS / "ecoli_strain_panel.csv",
                                  index=False)
    ok = [r for r in rows if r and "fisher_odds_ratio" in r]
    sig = [r for r in ok if r["fisher_p"] < 0.05]
    genes = {}
    for r in ok:
        for g in r["overrescued_genes"]:
            genes[g] = genes.get(g, 0) + 1
    summary = {"models_audited": len(ok), "models_significant": len(sig),
               "median_odds_ratio": float(pd.Series(
                   [r["fisher_odds_ratio"] for r in ok]).median()),
               "recurrent_genes": {g: n for g, n in
                                   sorted(genes.items(),
                                          key=lambda x: -x[1]) if n >= 3},
               "skipped": [r for r in rows if r and "fisher_odds_ratio"
                           not in r]}
    with open(RESULTS / "ecoli_strain_panel_summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)
    print(json.dumps(summary, indent=2)[:1500])


if __name__ == "__main__":
    main()
