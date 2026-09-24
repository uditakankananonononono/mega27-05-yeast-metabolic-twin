"""OMA ortholog-group recurrence of over-rescued genes across species.

If isozyme over-rescue is one conserved failure mode rather than six
model-specific artifacts, the same OMA ortholog groups should be hit in
multiple species. Each over-rescued gene's OMA group is fetched from the
OMA Browser REST API (yeast: systematic names; E. coli: UniProt
accessions from the committed uniprot_entries table; Mtb/B. subtilis:
locus tags; S. aureus: SA_RS via GFF old-tag bridge; Salmonella: STM).
Groups with over-rescued members in >=2 species are recurrences.
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

SETS = [
    ("yeast", "overrescue_audit.json"),
    ("ecoli", "overrescue_audit_ecoli.json"),
    ("bsubtilis", "overrescue_audit_bsubtilis.json"),
    ("mtb", "overrescue_audit_mtb.json"),
    ("saureus", "overrescue_audit_saureus_deg1002.json"),
    ("styphimurium", "overrescue_audit_styphimurium_deg1011.json"),
]


def oma_group(query: str):
    try:
        req = urllib.request.Request(
            f"https://omabrowser.org/api/protein/{query}/",
            headers={'User-Agent': 'Mozilla/5.0 (X11; Linux x86_64) '
                     'AppleWebKit/537.36 (KHTML, like Gecko) '
                     'Chrome/120.0 Safari/537.36',
                     'Accept': 'application/json'})
        with urllib.request.urlopen(req, timeout=30) as r:
            d = json.load(r)
        return d.get("oma_group"), d.get("omaid")
    except Exception:
        return None, None


def saureus_queries(genes):
    """SA_RS -> old SA tag via GFF (OMA indexes N315 as SA tags)."""
    m = {}
    for line in open(DATA / "gff_saureus_n315.gff"):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "gene":
            continue
        attrs = dict(re.findall(r"([^=;]+)=([^;]*)", f[8]))
        cur = attrs.get("locus_tag")
        for o in (attrs.get("old_locus_tag") or "").replace(
                "%2C", ",").split(","):
            if o:
                m.setdefault(cur, o)
    return [(g, m.get(g, g)) for g in genes]


def main():
    uni = pd.read_csv(RESULTS / "uniprot_entries_ecoli.csv")
    b2acc = dict(zip(uni.gene, uni.accession))
    rows = []
    for tag, audit in SETS:
        genes = json.load(open(RESULTS / audit))["overrescued_genes"]
        if tag == "ecoli":
            queries = [(g, b2acc.get(g, g)) for g in genes]
        elif tag == "saureus":
            queries = saureus_queries(genes)
        else:
            queries = [(g, g) for g in genes]
        for g, q in queries:
            grp, omaid = oma_group(q)
            rows.append({"species": tag, "gene": g, "query": q,
                         "oma_group": grp, "omaid": omaid})
            time.sleep(0.12)
    df = pd.DataFrame(rows)
    df.to_csv(RESULTS / "oma_groups_overrescued.csv", index=False)
    hit = df.dropna(subset=["oma_group"])
    rec = (hit.groupby("oma_group").species.nunique())
    rec = rec[rec >= 2].sort_values(ascending=False)
    detail = {int(g): sorted(hit[hit.oma_group == g]
                             .apply(lambda r: f"{r.species}:{r.gene}",
                                    axis=1)) for g in rec.index}
    out = {"n_genes": int(len(df)),
           "n_mapped": int(len(hit)),
           "n_groups_recurrent": int(len(rec)),
           "max_species_per_group": int(rec.max()) if len(rec) else 0,
           "recurrent_groups": detail}
    with open(RESULTS / "oma_recurrence.json", "w") as fh:
        json.dump(out, fh, indent=2)
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
