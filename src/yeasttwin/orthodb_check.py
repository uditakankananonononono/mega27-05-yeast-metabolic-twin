"""OrthoDB replication lane: paralog clustering on a second resource.

OMA groups place each over-rescued backup and its essential partner in
ortholog groups; OrthoDB at Eukaryota level additionally clusters
in-paralogs. Hypothesis: backup and partner fall in the same OrthoDB
Eukaryota-level group far above background, replicating the InterPro
domain-sharing result at the sequence-cluster level. Fisher exact test.
"""
import csv
import json
import time
import urllib.parse
import urllib.request
from itertools import combinations
from pathlib import Path

from scipy.stats import fisher_exact

ROOT = Path(__file__).resolve().parent.parent.parent
UA = {'User-Agent': 'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 '
                    '(KHTML, like Gecko) Chrome/120.0 Safari/537.36',
      'Accept': 'application/json'}


def orthodb_group(gene):
    q = urllib.parse.quote(gene)
    url = f"https://data.orthodb.org/current/search?query={q}&level=2759&limit=1"
    req = urllib.request.Request(url, headers=UA)
    try:
        with urllib.request.urlopen(req, timeout=30) as r:
            d = json.load(r)
        data = d.get("data") or []
        return data[0] if data else None
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
    grp = {}
    for g in genes:
        grp[g] = orthodb_group(g)
        time.sleep(0.2)
        print(g, grp[g], flush=True)
    pairset = {frozenset(p) for p in pairs}

    def same(a, b):
        return bool(grp[a] and grp[b] and grp[a] == grp[b])

    p_yes = sum(1 for a, b in pairs if same(a, b))
    bg = [pr for pr in combinations(genes, 2) if frozenset(pr) not in pairset]
    b_yes = sum(1 for a, b in bg if same(a, b))
    table = [[p_yes, len(pairs) - p_yes], [b_yes, len(bg) - b_yes]]
    orr, pval = fisher_exact(table)
    out_rows = [{"backup": a, "partner": b, "backup_group": grp[a] or "",
                 "partner_group": grp[b] or "", "same_group": same(a, b)}
                for a, b in pairs]
    with open(ROOT / "results" / "orthodb_pair_groups.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out_rows[0]))
        w.writeheader()
        w.writerows(out_rows)
    summary = {"n_genes": len(genes),
               "n_mapped": sum(1 for g in genes if grp[g]),
               "n_pairs": len(pairs), "pairs_same_group": p_yes,
               "n_background_pairs": len(bg), "background_same_group": b_yes,
               "fisher_oddsratio": orr, "fisher_p": pval}
    (ROOT / "results" / "orthodb_same_group.json").write_text(
        json.dumps(summary, indent=1))
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
