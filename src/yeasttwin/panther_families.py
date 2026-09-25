"""PANTHER classification lane for the over-rescued set.

Classifies the 19 over-rescued essential genes into PANTHER protein
families/subfamilies and PANTHER-GO-slim molecular functions via the
PANTHER web services (geneinfo), complementing the InterPro/Pfam domain
lane with family-level curation. Committed: results/panther_families.csv
(+ panther_summary.json).
"""
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
UA = {"User-Agent": "mega27-research/1.0"}
URL = "https://www.pantherdb.org/services/oai/pantherdb/geneinfo"


def geneinfo(orf):
    body = urllib.parse.urlencode(
        {"geneInputList": orf, "organism": "559292"}).encode()
    req = urllib.request.Request(URL, data=body, headers=UA)
    with urllib.request.urlopen(req, timeout=30) as r:
        d = json.load(r)
    g = d.get("search", {}).get("mapped_genes", {}).get("gene", None)
    if not g or not isinstance(g, dict):
        return {}
    fam, slim = "", ""
    adt = g.get("annotation_type_list", {}).get("annotation_data_type", [])
    if isinstance(adt, dict):
        adt = [adt]
    for block in adt:
        ann = block.get("annotation_list", {}).get("annotation", {})
        if isinstance(ann, dict):
            ann = [ann]
        names = ";".join(a.get("name", "") for a in ann)
        if block.get("content") == "ANNOT_TYPE_ID_PANTHER_GO_SLIM_MF":
            slim = names
    return {"family_id": g.get("family_id", ""), "subfamily_id": g.get("sf_id", ""),
            "go_slim_mf": slim}


def main():
    genes = [l.strip() for l in (ROOT / "results" / "overrescued_genes.csv")
             .read_text().splitlines() if l.strip()][1:]
    import csv
    rows = []
    for g in genes:
        try:
            info = geneinfo(g)
        except Exception as e:
            info = {"family_id": "", "subfamily_id": "", "go_slim_mf": f"ERROR {e}"}
        rows.append({"orf": g, **info})
        print(g, info.get("family_id"), flush=True)
        time.sleep(0.3)
    with open(ROOT / "results" / "panther_families.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    (ROOT / "results" / "panther_summary.json").write_text(json.dumps({
        "n_queried": len(rows),
        "n_classified": sum(1 for r in rows if r["family_id"]),
        "panther_version": 19}, indent=2))
    print(sum(1 for r in rows if r["family_id"]), "/", len(rows))


if __name__ == "__main__":
    main()
