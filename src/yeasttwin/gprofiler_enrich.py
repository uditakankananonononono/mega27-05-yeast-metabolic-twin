"""Functional enrichment of over-rescued gene sets (g:Profiler g:GOSt).

Characterizes WHICH functions get over-rescued in each species, and
contrasts yeast over-rescued essentials against FBA-caught essentials
(is over-rescue specific to isozyme-rich central carbon metabolism?).
One POST per gene set to the g:Profiler API; per-set significant terms
(Bonferroni-corrected by g:GOSt) written to results.
"""
from __future__ import annotations

import json
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results"

SETS = [
    ("yeast", "scerevisiae", "overrescue_audit.json"),
    ("ecoli", "ecoli", "overrescue_audit_ecoli.json"),
    ("bsubtilis", "bsubtilis", "overrescue_audit_bsubtilis.json"),
    ("mtb", "mtuberculosis", "overrescue_audit_mtb.json"),
]


def gost(query: list[str], organism: str) -> list[dict]:
    body = {"organism": organism, "query": query,
            "sources": ["GO:BP", "GO:MF", "KEGG", "REAC"],
            "no_evidences": True}
    req = urllib.request.Request(
        "https://biit.cs.ut.ee/gprofiler/api/gost/profile/",
        data=json.dumps(body).encode(),
        headers={"Content-Type": "application/json",
                 "User-Agent": "yeasttwin/1.0"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.load(r).get("result", [])


def main():
    summary = {}
    for tag, org, audit in SETS:
        genes = json.load(open(RESULTS / audit))["overrescued_genes"]
        try:
            res = gost(genes, org)
        except Exception as e:
            summary[tag] = {"error": str(e)[:80], "n_genes": len(genes)}
            print(tag, "ERROR", str(e)[:80])
            continue
        rows = [{"source": x["source"], "term": x["native"],
                 "name": x["name"], "p": x["p_value"],
                 "intersections": x.get("intersection_size")}
                for x in res]
        pd.DataFrame(rows).to_csv(
            RESULTS / f"gprofiler_overrescued_{tag}.csv", index=False)
        summary[tag] = {"organism": org, "n_genes": len(genes),
                        "n_significant_terms": len(res),
                        "top": rows[:5]}
        print(tag, len(genes), "genes ->", len(res), "terms;",
              rows[0]["name"] if rows else "none")
    # yeast contrast: caught essentials
    caught = json.loads((RESULTS / "overrescue_audit.json").read_text())
    n_caught = caught["n_caught"]
    # caught gene list is not stored; derive from labels + fba results
    from .labels import load_labels
    ko = pd.read_json(RESULTS / "fba_single_ko.json")["ratio"]
    ko = pd.Series(ko)
    labels = load_labels()
    ess = labels[labels].index
    pred_viable = ko.reindex(ess) > 1e-6
    caught_genes = sorted(ess[~pred_viable.fillna(False)])
    res = gost(caught_genes, "scerevisiae")
    pd.DataFrame([{"source": x["source"], "term": x["native"],
                   "name": x["name"], "p": x["p_value"]} for x in res]
                 ).to_csv(RESULTS / "gprofiler_caught_yeast.csv", index=False)
    summary["yeast_caught_contrast"] = {"n_genes": len(caught_genes),
                                        "n_significant_terms": len(res)}
    with open(RESULTS / "gprofiler_summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)
    print("yeast caught:", len(caught_genes), "->", len(res), "terms")


if __name__ == "__main__":
    main()
