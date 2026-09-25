"""PRIDE Archive lane: accession-level provenance for per-study abundances.

The PaxDb per-study abundance files used for the backup-expression
mechanism tests name their source studies (Ghaemmaghami, Newman, GPM, ...).
This lane queries the PRIDE Archive API for each yeast study token and
records matching PXD accessions, establishing accession-level proteomics
provenance for the abundance evidence. Committed:
results/pride_yeast_projects.csv (+ pride_provenance_summary.json).
"""
import json
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
UA = {"User-Agent": "mega27-research/1.0", "Accept": "application/json"}


def search(term):
    q = urllib.parse.quote(term)
    url = ("https://www.ebi.ac.uk/pride/ws/archive/v3/projects?"
           f"keyword={q}&pageSize=10&page=0")
    req = urllib.request.Request(url, headers=UA)
    try:
        with urllib.request.urlopen(req, timeout=30) as r:
            d = json.load(r)
    except Exception:
        return []
    return [{"pxd": p.get("accession", ""), "title": p.get("title", ""),
             "submissionDate": p.get("submissionDate", ""),
             "numAssays": p.get("numAssays", "")} for p in (d or [])]


def main():
    files = sorted((ROOT / "data" / "raw" / "paxdb_perstudy").glob("4932-*.txt"))
    tokens = []
    for f in files:
        t = re.sub(r"^4932-|\.txt$", "", f.name)
        t = re.sub(r"_(et_al_yeast_data.*|Saccharomyces.*)$", "", t)
        tokens.append((f.name, t))
    rows, n_with_pxd = [], 0
    for fname, tok in tokens:
        hits = search(tok) if re.search(r"[A-Za-z]{4,}", tok) else []
        rows.append({"paxdb_file": fname, "study_token": tok,
                     "n_pride_hits": len(hits),
                     "pxd_accessions": ";".join(h["pxd"] for h in hits),
                     "pxd_titles": " || ".join(h["title"] for h in hits)[:300]})
        n_with_pxd += bool(hits)
        print(tok, len(hits), flush=True)
        time.sleep(0.2)
    import csv
    with open(ROOT / "results" / "pride_yeast_projects.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    pxds = sorted({p for r in rows for p in r["pxd_accessions"].split(";") if p})
    (ROOT / "results" / "pride_provenance_summary.json").write_text(json.dumps({
        "n_paxdb_yeast_study_files": len(rows),
        "n_tokens_with_pride_hits": n_with_pxd,
        "n_distinct_pxd": len(pxds), "pxd_accessions": pxds}, indent=2))
    print(n_with_pxd, "/", len(rows), "tokens hit;", len(pxds), "PXD")


if __name__ == "__main__":
    main()
