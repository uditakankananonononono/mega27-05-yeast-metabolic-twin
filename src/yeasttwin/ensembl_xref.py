"""Ensembl REST lane: independent cross-reference of the over-rescued ORFs.

Verifies the 19 over-rescued ORF identifiers through Ensembl REST
(lookup/id), recording stable ID, display name, biotype and coordinates -
an identifier-resolution route independent of the SGD downloads and NCBI
lanes already used. Committed: results/ensembl_xref_overrescued.csv.
"""
import csv
import json
import time
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
UA = {"User-Agent": "mega27-research/1.0", "Content-Type": "application/json"}


def lookup(orf):
    url = f"https://rest.ensembl.org/lookup/id/{orf}?content-type=application/json"
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=30) as r:
        return json.load(r)


def main():
    genes = [l.strip() for l in (ROOT / "results" / "overrescued_genes.csv")
             .read_text().splitlines() if l.strip()][1:]
    rows = []
    for g in genes:
        try:
            d = lookup(g)
            rows.append({"orf": g, "ensembl_id": d.get("id", ""),
                         "display_name": d.get("display_name", ""),
                         "biotype": d.get("biotype", ""),
                         "seq_region": d.get("seq_region_name", ""),
                         "start": d.get("start", ""), "end": d.get("end", ""),
                         "source": d.get("source", "")})
        except Exception as e:  # recorded honest miss
            rows.append({"orf": g, "ensembl_id": "", "display_name": "",
                         "biotype": f"ERROR {e}", "seq_region": "",
                         "start": "", "end": "", "source": ""})
        print(g, rows[-1]["display_name"], flush=True)
        time.sleep(0.12)
    with open(ROOT / "results" / "ensembl_xref_overrescued.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    n_ok = sum(1 for r in rows if r["biotype"] == "protein_coding")
    (ROOT / "results" / "ensembl_xref_summary.json").write_text(json.dumps(
        {"n_queried": len(rows), "n_resolved_protein_coding": n_ok}, indent=2))
    print(n_ok, "/", len(rows))


if __name__ == "__main__":
    main()
