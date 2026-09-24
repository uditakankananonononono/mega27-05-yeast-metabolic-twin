"""UniProt annotation provenance for over-rescue pairs.

Fetches the curated UniProt entry for every over-rescued essential gene
and its model backup(s), per gene, from the UniProt REST API (yeast
organism 559292; E. coli K-12 organism 83333). Tests annotation transfer
at the curated-record level: (1) fraction of pairs where both members
carry only family-level SIMILARITY comments (evidence codes ECO:0000305,
sequence-similarity propagation); (2) asymmetry of CATALYTIC ACTIVITY
curation between essential and backup; (3) explicit non-complementation
language in backup records ("cannot", "does not", "no detectable",
"unable") - curated evidence the backup does not cover in vivo.
"""
from __future__ import annotations

import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

import cobra
import pandas as pd

from .paralogs import isozyme_partners

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data" / "raw"
RESULTS = ROOT / "results"
FIELDS = ("accession,id,gene_names,protein_name,cc_function,"
          "cc_similarity,cc_catalytic_activity")
NEG = ("cannot", "does not", "no detectable", "unable", "not essential",
       "dispensable", "not required")


def fetch_entry(gene: str, organism: int) -> dict:
    q = urllib.parse.quote(f"({gene}) AND (organism_id:{organism})")
    url = ("https://rest.uniprot.org/uniprotkb/search?query=" + q +
           "&fields=" + urllib.parse.quote(FIELDS) + "&format=tsv&size=1")
    with urllib.request.urlopen(url, timeout=30) as r:
        txt = r.read().decode()
    lines = txt.strip().split("\n")
    if len(lines) < 2:
        return {"gene": gene, "organism": organism}
    hdr, val = lines[0].split("\t"), lines[1].split("\t")
    rec = dict(zip(hdr, val))
    return {"gene": gene, "organism": organism,
            "accession": rec.get("Entry"), "entry_name": rec.get("Entry Name"),
            "protein": rec.get("Protein names", ""),
            "function": rec.get("Function [CC]", ""),
            "similarity": rec.get("Sequence similarities", ""),
            "catalytic": rec.get("Catalytic activity", "")}


def side(pairs: pd.DataFrame, organism: int, tag: str):
    genes = sorted(set(pairs.gene) | set(pairs.backup))
    recs = []
    for g in genes:
        try:
            recs.append(fetch_entry(g, organism))
        except Exception as e:
            recs.append({"gene": g, "organism": organism, "error": str(e)[:60]})
        time.sleep(0.15)
    df = pd.DataFrame(recs)
    df.to_csv(RESULTS / f"uniprot_entries_{tag}.csv", index=False)
    rows = []
    for _, r in pairs.iterrows():
        e = df[df.gene == r["gene"]].iloc[0]
        b = df[df.gene == r["backup"]].iloc[0]
        rows.append({
            "gene": r["gene"], "backup": r["backup"],
            "both_family_similarity": bool("SIMILARIT" in str(e.similarity)
                                           and "SIMILARIT" in str(b.similarity)),
            "essential_has_catalytic": bool(str(e.catalytic).strip()),
            "backup_has_catalytic": bool(str(b.catalytic).strip()),
            "backup_negative_language": any(k in (str(b.function)
                                                + str(b.catalytic)).lower()
                                            for k in NEG),
        })
    out = pd.DataFrame(rows)
    out.to_csv(RESULTS / f"uniprot_pairs_{tag}.csv", index=False)
    n = len(out)
    return {"tag": tag, "organism": organism, "n_pairs": n,
            "n_entries_fetched": int(len(df)),
            "both_family_similarity_frac": float(out.both_family_similarity.mean()),
            "essential_catalytic_frac": float(out.essential_has_catalytic.mean()),
            "backup_catalytic_frac": float(out.backup_has_catalytic.mean()),
            "backup_negative_language_frac": float(out.backup_negative_language.mean())}


def main():
    par = pd.read_csv(RESULTS / "overrescued_paralog_identity.csv")
    ypairs = pd.DataFrame(
        [{"gene": r["gene"], "backup": p}
         for _, r in par.iterrows()
         for p in str(r["partners"]).split(";") if p and p != "nan"])
    yres = side(ypairs, 559292, "yeast")
    print(json.dumps(yres, indent=2))

    emodel = cobra.io.load_json_model(str(DATA / "iML1515.json"))
    egenes = json.load(open(RESULTS / "overrescue_audit_ecoli.json"))[
        "overrescued_genes"]
    epairs = pd.DataFrame(
        [{"gene": g, "backup": b}
         for g in egenes for b in sorted(isozyme_partners(emodel, g))])
    eres = side(epairs, 83333, "ecoli")
    print(json.dumps(eres, indent=2))
    with open(RESULTS / "uniprot_pairs.json", "w") as fh:
        json.dump({"yeast": yres, "ecoli": eres}, fh, indent=2)


if __name__ == "__main__":
    main()
