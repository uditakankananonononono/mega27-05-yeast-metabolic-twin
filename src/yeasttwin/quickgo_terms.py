"""QuickGO (EBI GO annotation API) lane for the over-rescued set.

Fetches molecular-function GO annotations per UniProt accession for the 19
over-rescued essential genes and their backup partners, as annotation
evidence independent of SGD downloads and g:Profiler. Committed table:
results/quickgo_overrescued.csv (+ quickgo_summary.json).
"""
import csv
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent.parent
UA = {"User-Agent": "mega27-research/1.0", "Accept": "application/json"}


def go_mf(acc):
    q = urllib.parse.quote(f"geneProductId={acc}")
    url = ("https://www.ebi.ac.uk/QuickGO/services/annotation/search?"
           f"{q}&aspect=molecular_function&limit=100")
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=30) as r:
        d = json.load(r)
    out = {}
    for a in d.get("results", []):
        out[a["goId"]] = a.get("goName") or a["goId"]
    return out


def resolve_names(ids):
    out = {}
    for i in range(0, len(ids), 400):
        chunk = ids[i:i + 400]
        url = ("https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/"
               + ",".join(chunk))
        req = urllib.request.Request(url, headers=UA)
        with urllib.request.urlopen(req, timeout=60) as r:
            d = json.load(r)
        for t in d.get("results", []):
            out[t["id"]] = t.get("name") or t["id"]
    return out


def main():
    genes = [l.strip() for l in (ROOT / "results" / "overrescued_genes.csv")
             .read_text().splitlines() if l.strip()][1:]
    uni = pd.read_csv(ROOT / "results" / "uniprot_entries_yeast.csv")
    g2a = dict(zip(uni.gene, uni.accession))
    pairs = pd.read_csv(ROOT / "results" / "interpro_pair_domains.csv")
    b2a = dict(zip(pairs.backup, pairs.partner_acc))
    rows = []
    for g in genes:
        for role, acc in (("essential", g2a.get(g, "")),
                          ("backup", b2a.get(g, ""))):
            if not acc:
                continue
            terms = go_mf(acc)
            rows.append({"gene": g, "role": role, "accession": acc,
                         "n_mf_terms": len(terms),
                         "go_ids": ";".join(sorted(terms)),
                         "go_names": ""})
            print(g, role, acc, len(terms), flush=True)
            time.sleep(0.15)
    df = pd.DataFrame(rows)
    names = resolve_names(sorted({g for t in df.go_ids for g in t.split(";") if g}))
    df["go_names"] = df.go_ids.map(lambda s: ";".join(sorted({names.get(g, g) for g in s.split(";") if g})))
    df.to_csv(ROOT / "results" / "quickgo_overrescued.csv", index=False)
    cat = df.go_names.str.contains("catalytic|activity", case=False, na=False)
    summary = {
        "n_genes": int(len(genes)),
        "n_records": int(len(df)),
        "essential_median_mf_terms": float(df[df.role == "essential"].n_mf_terms.median()),
        "backup_median_mf_terms": float(df[df.role == "backup"].n_mf_terms.median()),
        "frac_with_activity_term": float(cat.mean()),
    }
    (ROOT / "results" / "quickgo_summary.json").write_text(json.dumps(summary, indent=2))
    print(summary)


if __name__ == "__main__":
    main()
