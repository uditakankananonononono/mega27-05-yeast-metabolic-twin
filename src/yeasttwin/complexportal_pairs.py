"""EBI Complex Portal lane: co-complex membership of over-rescue pairs.

Queries the Complex Portal API per UniProt accession for each
essential/backup pair and tests whether the two proteins share a curated
macromolecular complex - a physical-mechanism lane orthogonal to the
genetic-interaction (BioGRID) and domain (InterPro) evidence.
Committed: results/complexportal_pairs.csv (+ complexportal_summary.json).
"""
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent.parent
UA = {"User-Agent": "mega27-research/1.0", "Accept": "application/json"}


def complexes(acc):
    url = ("https://www.ebi.ac.uk/complexportal/api/search/"
           f"{urllib.parse.quote(acc)}?format=json&first=0&number=25")
    req = urllib.request.Request(url, headers=UA)
    try:
        with urllib.request.urlopen(req, timeout=30) as r:
            d = json.load(r)
    except Exception:
        return {}
    out = {}
    for hit in d.get("hits", []):
        j = hit.get("_source", {})
        ac = j.get("complex_ac") or hit.get("_id", "")
        out[ac] = j.get("complex_name", "")
    return out


def main():
    pairs = pd.read_csv(ROOT / "results" / "interpro_pair_domains.csv")
    rows = []
    for _, r in pairs.iterrows():
        ca, cb = complexes(r.backup_acc), complexes(r.partner_acc)
        shared = sorted(set(ca) & set(cb))
        rows.append({"essential": r.backup, "backup": r.partner,
                     "essential_acc": r.backup_acc, "backup_acc": r.partner_acc,
                     "n_complexes_essential": len(ca), "n_complexes_backup": len(cb),
                     "shared_complexes": ";".join(shared),
                     "shared_complex_names": ";".join(ca[s] for s in shared),
                     "co_complex": bool(shared)})
        print(r.backup, r.partner, bool(shared), flush=True)
        time.sleep(0.15)
    df = pd.DataFrame(rows)
    df.to_csv(ROOT / "results" / "complexportal_pairs.csv", index=False)
    (ROOT / "results" / "complexportal_summary.json").write_text(json.dumps({
        "n_pairs": int(len(df)),
        "n_co_complex": int(df.co_complex.sum()),
        "note": "honest null retained when no shared curated complex"}, indent=2))
    print(df.co_complex.sum(), "/", len(df))


if __name__ == "__main__":
    main()
