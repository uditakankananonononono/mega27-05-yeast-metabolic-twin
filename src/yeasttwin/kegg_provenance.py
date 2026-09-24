"""KEGG provenance check: do over-rescued essentials and their model
backups share KEGG ortholog (KO) annotations?

Shared KO = both genes received the same functional annotation, the
likely route by which the model's GPR acquired the false isozyme clause.
Uses the KEGG REST API; results cached to data/raw/kegg_ko_cache.tsv.
"""
from __future__ import annotations

import json
import re
import time
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data" / "raw"
RESULTS = ROOT / "results"
CACHE = DATA / "kegg_ko_cache.tsv"


def fetch_ko(genes: list[str]) -> dict[str, str]:
    """ORF -> KO id (or '' if none). Cached; 10 genes per API call."""
    cache = {}
    if CACHE.exists():
        cache = pd.read_csv(CACHE, sep="\t", dtype=str).set_index("orf")["ko"].to_dict()
    todo = [g for g in genes if g not in cache]
    for i in range(0, len(todo), 10):
        batch = todo[i:i + 10]
        url = "https://rest.kegg.jp/get/" + "+".join("sce:" + g for g in batch)
        txt = urllib.request.urlopen(url, timeout=30).read().decode()
        for entry in txt.split("///"):
            m = re.search(r"ENTRY\s+(Y[A-Z0-9]{5,6})\b", entry)
            k = re.search(r"ORTHOLOGY\s+(K\d{5})", entry)
            if m:
                cache[m.group(1)] = k.group(1) if k else ""
        time.sleep(0.4)
        pd.DataFrame({"orf": list(cache), "ko": list(cache.values())}
                     ).to_csv(CACHE, sep="\t", index=False)
    return {g: cache.get(g, "") for g in genes}


def main():
    par = pd.read_csv(RESULTS / "overrescued_paralog_identity.csv")
    genes = set(par.gene)
    for ps in par.partners.dropna():
        genes.update(str(ps).split(";"))
    ko = fetch_ko(sorted(genes))
    rows = []
    for _, r in par.iterrows():
        partners = [p for p in str(r["partners"]).split(";") if p and p != "nan"]
        same = [p for p in partners if ko.get(p) and ko.get(p) == ko.get(r["gene"])]
        rows.append({"gene": r["gene"], "ko": ko.get(r["gene"], ""),
                     "n_partners": len(partners),
                     "n_same_ko": len(same),
                     "same_ko_partners": ";".join(same)})
    df = pd.DataFrame(rows)
    df.to_csv(RESULTS / "kegg_same_ko_backups.csv", index=False)
    tot_p = int(df.n_partners.sum())
    same_p = int(df.n_same_ko.sum())
    out = {"n_genes": len(df), "n_backup_pairs": tot_p,
           "n_same_ko_pairs": same_p,
           "frac_same_ko": same_p / max(tot_p, 1)}
    with open(RESULTS / "kegg_same_ko_backups.json", "w") as fh:
        json.dump(out, fh, indent=2)
    print(json.dumps(out, indent=2))
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
