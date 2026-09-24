"""InterPro/Pfam domain-sharing lane.

For each over-rescued backup and its essential partner (same KEGG KO), fetch
the UniProt accession and the protein's Pfam domains from the InterPro API.
Hypothesis (annotation-transfer mechanism): backups share Pfam domains with
their partners far above the background rate among the other gene pairs in
the same protein set. Fisher exact test on shared-domain yes/no.
"""
import csv
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

from scipy.stats import fisher_exact
from itertools import combinations

ROOT = Path(__file__).resolve().parent.parent.parent
UA = {'User-Agent': 'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 '
                    '(KHTML, like Gecko) Chrome/120.0 Safari/537.36',
      'Accept': 'application/json'}


def get_json(url):
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=30) as r:
        return json.load(r)


def uniprot_acc(gene):
    q = urllib.parse.quote(f"{gene} AND organism_id:559292")
    url = ("https://rest.uniprot.org/uniprotkb/search?query=" + q +
           "&fields=accession,gene_names&format=tsv&size=1")
    req = urllib.request.Request(url, headers=UA)
    try:
        with urllib.request.urlopen(req, timeout=30) as r:
            lines = r.read().decode().strip().splitlines()
        return lines[1].split("\t")[0] if len(lines) > 1 else None
    except Exception:
        return None


def pfam_domains(acc):
    try:
        d = get_json("https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/"
                     f"uniprot/{acc}/?page_size=200")
        return sorted({e["metadata"]["accession"] for e in d.get("results", [])})
    except Exception:
        return []


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
    acc, dom = {}, {}
    for g in genes:
        acc[g] = uniprot_acc(g)
        time.sleep(0.15)
        dom[g] = pfam_domains(acc[g]) if acc[g] else []
        time.sleep(0.15)
        print(g, acc[g], len(dom[g]), flush=True)
    pairset = {frozenset(p) for p in pairs}

    def shared(a, b):
        return bool(set(dom[a]) & set(dom[b]))

    p_yes = sum(1 for a, b in pairs if shared(a, b))
    bg = [pr for pr in combinations(genes, 2) if frozenset(pr) not in pairset]
    b_yes = sum(1 for a, b in bg if shared(a, b))
    table = [[p_yes, len(pairs) - p_yes], [b_yes, len(bg) - b_yes]]
    orr, pval = fisher_exact(table)
    out_rows = [{"backup": a, "partner": b,
                 "backup_acc": acc[a], "partner_acc": acc[b],
                 "backup_domains": ";".join(dom[a]),
                 "partner_domains": ";".join(dom[b]),
                 "shared_domains": ";".join(sorted(set(dom[a]) & set(dom[b]))),
                 "shares_domain": shared(a, b)} for a, b in pairs]
    with open(ROOT / "results" / "interpro_pair_domains.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out_rows[0]))
        w.writeheader()
        w.writerows(out_rows)
    summary = {"n_genes": len(genes), "n_mapped": sum(1 for g in genes if acc[g]),
               "n_pairs": len(pairs), "pairs_sharing_domain": p_yes,
               "n_background_pairs": len(bg), "background_sharing_domain": b_yes,
               "fisher_oddsratio": orr, "fisher_p": pval}
    (ROOT / "results" / "interpro_domain_sharing.json").write_text(
        json.dumps(summary, indent=1))
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
