"""RCSB PDB structural coverage of over-rescue pairs.

For each yeast over-rescued essential gene and its backup(s): does an
experimental structure exist (via UniProt accession -> RCSB search API)?
A backup with a propagated catalytic annotation but NO experimental
structure is pure in-silico annotation transfer; structural coverage
asymmetry between pair members is recorded per pair.
"""
from __future__ import annotations

import json
import time
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results"


def rcsb_count(uniprot_acc: str) -> tuple[int, str]:
    q = {"query": {"type": "terminal", "service": "text", "parameters": {
        "attribute": "rcsb_polymer_entity_container_identifiers"
                     ".reference_sequence_identifiers.database_accession",
        "operator": "exact_match", "value": uniprot_acc}},
        "return_type": "entry",
        "request_options": {"paginate": {"start": 0, "rows": 3}}}
    req = urllib.request.Request(
        "https://search.rcsb.org/rcsbsearch/v2/query",
        data=json.dumps(q).encode(),
        headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=30) as r:
        d = json.load(r)
    ids = [x["identifier"] for x in d.get("result_set", [])]
    return d.get("total_count", 0), ";".join(ids)


def main():
    entries = pd.read_csv(RESULTS / "uniprot_entries_yeast.csv")
    pairs = pd.read_csv(RESULTS / "uniprot_pairs_yeast.csv")
    acc = dict(zip(entries.gene, entries.accession))
    genes = sorted(set(pairs.gene) | set(pairs.backup))
    cover = {}
    for g in genes:
        a = acc.get(g)
        if not a or pd.isna(a):
            cover[g] = (0, "")
            continue
        try:
            cover[g] = rcsb_count(a)
        except Exception:
            cover[g] = (0, "")
        time.sleep(0.2)
    rows = [{"gene": g, "uniprot": acc.get(g),
             "n_pdb": cover[g][0], "pdb_ids": cover[g][1]} for g in genes]
    df = pd.DataFrame(rows)
    df.to_csv(RESULTS / "pdb_structures_yeast_pairs.csv", index=False)
    ess = set(pairs.gene); bak = set(pairs.backup)
    out = {
        "n_genes": len(genes),
        "essential_with_structure": int(sum(cover[g][0] > 0 for g in ess)),
        "n_essential": len(ess),
        "backup_with_structure": int(sum(cover[g][0] > 0 for g in bak)),
        "n_backup": len(bak),
        "backups_without_structure": sorted(g for g in bak
                                            if cover[g][0] == 0),
    }
    with open(RESULTS / "pdb_pairs.json", "w") as fh:
        json.dump(out, fh, indent=2)
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
