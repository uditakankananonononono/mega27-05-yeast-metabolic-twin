"""EBI OLS4 / ChEBI ontology lane for design metabolites.

Looks up exact ChEBI terms through the EBI Ontology Lookup Service for the
metabolites of the sdh2/ach1 succinate design, recording the ChEBI ID,
label and description used in the design rationale.
"""
import csv
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
UA = {'User-Agent': 'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 '
                    '(KHTML, like Gecko) Chrome/120.0 Safari/537.36',
      'Accept': 'application/json'}
TERMS = ["succinate", "fumarate", "acetyl-CoA", "coenzyme A", "ubiquinone"]


def ols_exact(term):
    q = urllib.parse.quote(term)
    url = (f"https://www.ebi.ac.uk/ols4/api/search?q={q}"
           f"&ontology=chebi&exact=true&rows=1")
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=30) as r:
        d = json.load(r)
    docs = d.get("response", {}).get("docs", [])
    return docs[0] if docs else {}


def main():
    rows = []
    for t in TERMS:
        d = ols_exact(t)
        rows.append({"query": t, "chebi_id": d.get("obo_id", ""),
                     "label": d.get("label", ""),
                     "description": (d.get("description") or [""])[0]})
        print(t, rows[-1]["chebi_id"], flush=True)
        time.sleep(0.2)
    with open(ROOT / "results" / "ols_chebi_design_metabolites.csv",
              "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()
