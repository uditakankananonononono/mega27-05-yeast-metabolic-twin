"""IntAct (EBI PSICQUIC) lane: physical interaction evidence for pairs.

Counts curated physical interactions between each essential gene's UniProt
accession and its backup's accession via the IntAct PSICQUIC service
(format=count) - orthogonal to the BioGRID genetic-interaction lane, which
covers genetic rescue rather than physical binding.
Committed: results/intact_pairs.csv (+ intact_summary.json).
"""
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent.parent
UA = {"User-Agent": "mega27-research/1.0"}
BASE = ("https://www.ebi.ac.uk/Tools/webservices/psicquic/intact/"
        "webservices/current/search/query/")


def count_pair(a, b):
    q = urllib.parse.quote(f"{a} AND {b}")
    req = urllib.request.Request(BASE + q + "?format=count", headers=UA)
    try:
        with urllib.request.urlopen(req, timeout=30) as r:
            return int(r.read().decode().strip() or 0)
    except Exception:
        return -1


def main():
    pairs = pd.read_csv(ROOT / "results" / "interpro_pair_domains.csv")
    rows = []
    for _, r in pairs.iterrows():
        n = count_pair(r.backup_acc, r.partner_acc)
        rows.append({"essential": r.backup, "backup": r.partner,
                     "essential_acc": r.backup_acc, "backup_acc": r.partner_acc,
                     "intact_physical_interactions": n})
        print(r.backup, r.partner, n, flush=True)
        time.sleep(0.2)
    df = pd.DataFrame(rows)
    df.to_csv(ROOT / "results" / "intact_pairs.csv", index=False)
    ok = df[df.intact_physical_interactions >= 0]
    (ROOT / "results" / "intact_summary.json").write_text(json.dumps({
        "n_pairs": int(len(df)),
        "n_with_physical_evidence": int((ok.intact_physical_interactions > 0).sum()),
        "n_query_errors": int((df.intact_physical_interactions < 0).sum()),
        "note": "physical binding is not expected for rescue-by-isozyme; "
                "honest low counts retained"}, indent=2))
    print(ok.intact_physical_interactions.sum())


if __name__ == "__main__":
    main()
