"""RDKit lane: descriptor computation for the succinate-design metabolites.

Recomputes molecular descriptors (exact mass, TPSA, H-bond donors/acceptors,
rotatable bonds, formal charge, fraction Csp3) from PubChem isomeric SMILES
with RDKit, cross-checking the committed PubChem property records used in
the design's redox/charge rationale. Discrepancies are recorded honestly.
Committed: results/rdkit_design_metabolites.csv (+ rdkit_summary.json).
"""
import json
import time
import urllib.request
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

ROOT = Path(__file__).resolve().parent.parent.parent
UA = {"User-Agent": "mega27-research/1.0", "Accept": "application/json"}


def smiles_for(cid):
    url = (f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}"
           "/property/IsomericSMILES,ConnectivitySMILES/JSON")
    req = urllib.request.Request(url, headers=UA)
    d = None
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req, timeout=30) as r:
                d = json.load(r)
            break
        except Exception:
            time.sleep(3 * (attempt + 1))
    if d is None:
        raise RuntimeError(f"PubChem unreachable for CID {cid}")
    p = d["PropertyTable"]["Properties"][0]
    return p.get("IsomericSMILES") or p.get("ConnectivitySMILES")


def main():
    pub = pd.read_csv(ROOT / "results" / "pubchem_tca_compounds.csv")
    rows = []
    for _, r in pub.iterrows():
        smi = smiles_for(int(r.CID))
        m = Chem.MolFromSmiles(smi)
        rows.append({
            "query": r.query, "CID": int(r.CID), "smiles": smi,
            "rdkit_exact_mw": round(Descriptors.ExactMolWt(m), 3),
            "rdkit_tpsa": round(rdMolDescriptors.CalcTPSA(m), 2),
            "rdkit_hbd": rdMolDescriptors.CalcNumHBD(m),
            "rdkit_hba": rdMolDescriptors.CalcNumHBA(m),
            "rdkit_rotatable": rdMolDescriptors.CalcNumRotatableBonds(m),
            "rdkit_formal_charge": Chem.GetFormalCharge(m),
            "rdkit_fraction_csp3": round(rdMolDescriptors.CalcFractionCSP3(m), 3),
            "pubchem_mw": r.MolecularWeight, "pubchem_tpsa": r.TPSA,
            "tpsa_agree": abs(rdMolDescriptors.CalcTPSA(m) - float(r.TPSA)) < 0.51})
        print(r.query, smi, flush=True)
        time.sleep(0.2)
    df = pd.DataFrame(rows)
    df.to_csv(ROOT / "results" / "rdkit_design_metabolites.csv", index=False)
    (ROOT / "results" / "rdkit_summary.json").write_text(json.dumps({
        "n_compounds": int(len(df)),
        "n_tpsa_agree_pubchem": int(df.tpsa_agree.sum()),
        "dianion_set": df.loc[df.rdkit_formal_charge <= -2, "query"].tolist()}, indent=2))
    print(df[["query", "rdkit_formal_charge", "tpsa_agree"]])


if __name__ == "__main__":
    main()
