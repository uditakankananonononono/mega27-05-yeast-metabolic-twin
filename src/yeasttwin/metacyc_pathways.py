"""BioCyc/MetaCyc pathway co-membership lane.

For each over-rescued backup and its essential partner, fetch the yeast
BioCyc (YEAST, built on MetaCyc) pathway memberships via the BioCyc web
services apixml interface. Hypothesis: backups and partners co-occur in at
least one curated pathway far above the background rate among the other
gene pairs in the same set (the rescue is within-pathway, not distant).
Fisher exact test on co-membership yes/no.
"""
import csv
import json
import re
import time
import urllib.request
from itertools import combinations
from pathlib import Path

from scipy.stats import fisher_exact

ROOT = Path(__file__).resolve().parent.parent.parent
UA = {'User-Agent': 'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 '
                    '(KHTML, like Gecko) Chrome/120.0 Safari/537.36'}
FRAME_RE = re.compile(r"<Pathway ID='YEAST:([^']+)'")


def pathways_of_gene(orf):
    url = (f"https://websvc.biocyc.org/apixml?fn=pathways-of-gene"
           f"&id=YEAST:{orf}")
    req = urllib.request.Request(url, headers=UA)
    try:
        with urllib.request.urlopen(req, timeout=30) as r:
            return sorted(set(FRAME_RE.findall(r.read().decode("utf-8", "replace"))))
    except Exception:
        return None


def main():
    rows = list(csv.DictReader(
        open(ROOT / "results" / "kegg_same_ko_backups.csv")))
    pairs = set()
    for r in rows:
        for p in r["same_ko_partners"].split(";"):
            if p:
                pairs.add((r["gene"], p))
    pairs = sorted(pairs)
    genes = sorted({g for pr in pairs for g in pr})
    pw = {}
    for g in genes:
        pw[g] = pathways_of_gene(g)
        time.sleep(0.2)
        print(g, None if pw[g] is None else len(pw[g]), flush=True)
    pairset = {frozenset(p) for p in pairs}

    def co(a, b):
        return bool(pw[a] and pw[b] and set(pw[a]) & set(pw[b]))

    p_yes = sum(1 for a, b in pairs if co(a, b))
    bg = [pr for pr in combinations(genes, 2) if frozenset(pr) not in pairset]
    b_yes = sum(1 for a, b in bg if co(a, b))
    table = [[p_yes, len(pairs) - p_yes], [b_yes, len(bg) - b_yes]]
    orr, pval = fisher_exact(table)
    out_rows = [{"backup": a, "partner": b,
                 "backup_pathways": ";".join(pw[a] or []),
                 "partner_pathways": ";".join(pw[b] or []),
                 "shared_pathways": ";".join(sorted(set(pw[a] or []) & set(pw[b] or []))),
                 "co_member": co(a, b)} for a, b in pairs]
    with open(ROOT / "results" / "metacyc_pair_pathways.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out_rows[0]))
        w.writeheader()
        w.writerows(out_rows)
    summary = {"n_genes": len(genes),
               "n_mapped": sum(1 for g in genes if pw[g] is not None),
               "n_pairs": len(pairs), "pairs_co_member": p_yes,
               "n_background_pairs": len(bg), "background_co_member": b_yes,
               "fisher_oddsratio": orr, "fisher_p": pval}
    (ROOT / "results" / "metacyc_pathway_summary.json").write_text(
        json.dumps(summary, indent=1))
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
