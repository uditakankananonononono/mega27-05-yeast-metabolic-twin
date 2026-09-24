"""Rhea reaction-record lane for the sdh2/ach1 succinate design.

Retrieves the curated Rhea reactions for the design enzymes (SDH2:
succinate dehydrogenase ubiquinone EC 1.3.5.1; ACH1: acetyl-CoA hydrolase
EC 3.1.2.1) and the main competing acetyl-CoA sink (citrate synthase
EC 2.3.3.1), confirming the deletion targets and the obligatory
succinate-yielding route against an expert-curated reaction resource.
"""
import csv
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
UA = {'User-Agent': 'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 '
                    '(KHTML, like Gecko) Chrome/120.0 Safari/537.36'}
QUERIES = [("SDH2", "1.3.5.1", "design knockout: succinate dehydrogenase"),
           ("ACH1", "3.1.2.1", "design knockout: acetyl-CoA hydrolase"),
           ("CIT", "2.3.3.1", "competing acetyl-CoA sink: citrate synthase")]


def rhea_by_ec(ec, limit=5):
    q = urllib.parse.quote(f"ec:{ec}")
    url = (f"https://www.rhea-db.org/rhea?query={q}"
           f"&columns=rhea-id,equation,ec&format=tsv&limit={limit}")
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=30) as r:
        lines = r.read().decode().strip().splitlines()
    return [dict(zip(["rhea_id", "equation", "ec"],
                     ln.split("\t") + [""] * (3 - len(ln.split("\t")))))
            for ln in lines[1:]]


def main():
    rows = []
    for gene, ec, role in QUERIES:
        for rec in rhea_by_ec(ec):
            rec.update({"gene": gene, "role": role})
            rows.append(rec)
        print(gene, ec, len(rows), flush=True)
    with open(ROOT / "results" / "rhea_design_reactions.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["gene", "ec", "role", "rhea_id", "equation"])
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()
