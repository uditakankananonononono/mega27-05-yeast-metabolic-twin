"""Accession-level dataset ledger with honest source-study deduplication.

Scans the committed result/data files and emits one row per accession-level
dataset record actually used in the analyses (UniProt entries, PDB IDs,
InterPro domains, OMA/OrthoDB groups, PaxDb per-study files, RefSeq GFFs,
Expression Atlas experiments, Rhea/ChEBI/PubChem records, PRIDE projects).

Collapse rules (applied, recorded in the summary):
- exact duplicate (resource, accession) rows count once;
- integrated PaxDb whole-organism files are NOT counted when their per-study
  files are counted (re-deposit of the same studies at coarser granularity);
- per-species repeated resource releases (e.g. the four E. coli model
  generations) are distinct studies and count individually;
- literature records (Crossref DOIs) are NOT datasets and are excluded.

Committed: results/accession_ledger.csv (+ accession_ledger_summary.json).
"""
import csv
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent.parent
RES = ROOT / "results"
RAW = ROOT / "data" / "raw"


def add(rows, resource, accession, study, used_for, evidence):
    accession = str(accession).strip()
    if accession:
        rows.append({"resource": resource, "accession": accession,
                     "source_study_or_release": study, "used_for": used_for,
                     "evidence_file": evidence})


def main():
    rows = []
    uni_y = pd.read_csv(RES / "uniprot_entries_yeast.csv")
    for _, r in uni_y.iterrows():
        add(rows, "UniProtKB", r.accession, "UniProtKB/Swiss-Prot release 2026-09",
            "over-rescue pair annotation (yeast)", "results/uniprot_entries_yeast.csv")
    uni_e = pd.read_csv(RES / "uniprot_entries_ecoli.csv")
    for _, r in uni_e.iterrows():
        add(rows, "UniProtKB", r.accession, "UniProtKB/Swiss-Prot release 2026-09",
            "over-rescue pair annotation (E. coli)", "results/uniprot_entries_ecoli.csv")
    pdb = pd.read_csv(RES / "pdb_structures_yeast_pairs.csv")
    for _, r in pdb.iterrows():
        for pid in str(r.pdb_ids).split(";"):
            add(rows, "RCSB PDB", pid, "RCSB PDB release 2026-09",
                f"structural coverage of pair ({r.gene})",
                "results/pdb_structures_yeast_pairs.csv")
    ipr = pd.read_csv(RES / "interpro_pair_domains.csv")
    for _, r in ipr.iterrows():
        for dom in set(str(r.backup_domains).split(";") + str(r.partner_domains).split(";")):
            add(rows, "InterPro (Pfam)", dom, "InterPro release 2026-09",
                "domain-sharing test of backup pairs", "results/interpro_pair_domains.csv")
    oma = pd.read_csv(RES / "oma_groups_overrescued.csv")
    for _, r in oma.iterrows():
        add(rows, "OMA Browser", r.omaid, "OMA release 2026-09",
            f"ortholog-group recurrence ({r.species}:{r.gene})",
            "results/oma_groups_overrescued.csv")
    odb = pd.read_csv(RES / "orthodb_pair_groups.csv")
    for col in ["backup_group", "partner_group"]:
        for v in odb[col].unique():
            add(rows, "OrthoDB", v, "OrthoDB v12", "sequence-cluster paralogy check",
                "results/orthodb_pair_groups.csv")
    rhea = pd.read_csv(RES / "rhea_design_reactions.csv")
    for _, r in rhea.iterrows():
        for c in rhea.columns:
            if "rhea" in c.lower() or "id" == c.lower():
                add(rows, "Rhea", r[c], "Rhea release 2026-09", "design reaction curation",
                    "results/rhea_design_reactions.csv")
                break
    chebi = pd.read_csv(RES / "ols_chebi_design_metabolites.csv")
    for _, r in chebi.iterrows():
        add(rows, "ChEBI via EBI OLS4", r.chebi_id, "ChEBI release via OLS4",
            "design metabolite ontology", "results/ols_chebi_design_metabolites.csv")
    pub = pd.read_csv(RES / "pubchem_tca_compounds.csv")
    for _, r in pub.iterrows():
        add(rows, "PubChem", f"CID{r.CID}", "PubChem release 2026-09",
            "design redox/charge context", "results/pubchem_tca_compounds.csv")
    for f in sorted(RAW.glob("gff_*.gff")):
        acc = {"gff_bsubtilis_168.gff": "GCF_000009045.1",
               "gff_ecoli_mg1655.gff": "GCF_000005845.2",
               "gff_hpylori_26695.gff": "GCF_000008525.1",
               "gff_mtb_h37rv.gff": "GCF_000195955.2",
               "gff_salmonella_lt2.gff": "GCF_000006945.2",
               "gff_saureus_n315.gff": "GCF_000009645.1"}[f.name]
        add(rows, "NCBI RefSeq (Datasets API)", acc, "RefSeq annotation release",
            "DEG-to-model gene mapping", f"data/raw/{f.name}")
    for f in sorted((RAW / "paxdb_perstudy").glob("*.txt")):
        if "WHOLE_ORGANISM-integrated" in f.name:
            continue  # re-deposit of the per-study files at coarser granularity
        species = f.name.split("-")[0]
        add(rows, "PaxDb per-study proteomics", f.name.replace(".txt", ""),
            f"PaxDb v4.2 per-study dataset (taxid {species})",
            "per-study abundance replication", f"data/raw/paxdb_perstudy/{f.name}")
    for f in sorted((RAW / "gxa_yeast").glob("E-MTAB-*.tsv")):
        add(rows, "EBI Expression Atlas", f.name.split("-tpms")[0],
            "Expression Atlas baseline RNA-seq", "isozyme co-expression test",
            f"data/raw/gxa_yeast/{f.name}")
    pride = RES / "pride_yeast_projects.csv"
    if pride.exists():
        pr = pd.read_csv(pride).fillna("")
        for _, r in pr.iterrows():
            for pxd in str(r.pxd_accessions).split(";"):
                add(rows, "PRIDE Archive", pxd, "PRIDE Archive 2026-09",
                    f"provenance for PaxDb study {r.study_token}",
                    "results/pride_yeast_projects.csv")
    df = pd.DataFrame(rows).drop_duplicates(subset=["resource", "accession"])
    df.to_csv(RES / "accession_ledger.csv", index=False)
    by_res = df.groupby("resource").size().to_dict()
    summary = {
        "total_distinct_accession_level_records": int(len(df)),
        "by_resource": {k: int(v) for k, v in sorted(by_res.items())},
        "collapse_rules": [
            "duplicate (resource, accession) rows count once",
            "PaxDb whole-organism integrated files excluded as re-deposits of the counted per-study files",
            "distinct source studies count individually even within one resource",
            "Crossref/literature DOIs excluded (not datasets)"],
        "gate_120": bool(len(df) >= 120),
    }
    (RES / "accession_ledger_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
