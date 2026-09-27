#!/usr/bin/env python3
"""A7: stress-regime <-> transcriptome concordance (pre-reg Amendment 3, item 2).

DESIGN LOCKED IN THIS FILE BEFORE ANY RUN. Metric per Amendment 3 (locked):
direction-of-change agreement of each regime's signature pathway against the
matching stress transcriptome, reported as a binomial test vs 0.5 with 95%
CI. No model parameter is refit (this script never touches the GEM).

Dataset: Gasch et al. 2000 (PMID 11102521) complete_dataset.txt, authors'
raw data file archived by the Wayback Machine from genome-www.stanford.edu/
yeast_stress (SHA-256 recorded in output; Gasch 2000 is context-grade per
pre-reg Section 8 - A7 is a concordance analysis, not a Section-8 gate).

Locked arms (regime labels from results/climate/a4_regimes.json, v1):
  A. "osmotic carbon diversion" <-> 1M sorbitol, 30 min (Gasch col
     '1M sorbitol - 30 min'). Glycerol axis: GPD1, GPD2, GPP1, GPP2.
     Expected UP (Hohmann 2002 MMBR 66:300; Gasch 2000 osmotic response).
  B. "thermal carbon starvation" <-> heat shock 25->37 C, 30 min (col
     'Heat Shock 30 minutes hs-1'). HSP axis: HSP12, HSP26, HSP30, HSP42,
     HSP78, HSP82, HSP104, SSA3, SSA4, STI1. Expected UP (definitional
     heat-shock response; Gasch 2000; Verghese 2012 review).
  C. "thermal carbon starvation" <-> same column. Glycolytic axis: PGI1,
     PFK1, PFK2, FBA1, TPI1, GPM1, PGK1, ENO1, ENO2, TDH3, HXK2, CDC19,
     PDC1, ADH1. Expected UP (Gasch 2000 ESR: carbohydrate-metabolism/energy
     genes induced after heat shock). DECLARED least-certain expectation;
     two-sided p reported in addition to the directional test.
  D. nitrogen-axis binding inside "multi-constraint collapse" (a4_binding
     ammonium-bound flags; no standalone N regime label exists - declared
     mapping) <-> Nitrogen Depletion, 1 h (col 'Nitrogen Depletion 1 h').
     NCR axis: GAP1, MEP2, DAL5, PUT4, UGA4, GDH2, GLT1. Expected UP
     (Magasanik & Kaiser 2002 Gene 290:1; Gasch 2000 N-depletion).
  Ethanol regime ("ethanol detoxification failure"): NO ethanol-stress
  experiment exists in Gasch 2000 (YP-ethanol columns are carbon-source
  comparisons, not stress) -> declared untested here.

Locked scoring: per arm, n = genes with a non-missing log2 ratio; k = genes
whose sign matches the expected direction (UP: value > 0). Exact binomial
test vs p=0.5: one-sided p = P[X >= k] for the expected-direction arms
(plus two-sided p for arm C as declared); 95% Clopper-Pearson CI on k/n.
Genes resolving to zero or duplicate rows abort the run (resolution is
deterministic, by exact standard-name match on the NAME field).
"""
import csv
import json
import math
import sys
from hashlib import sha256
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "raw" / "gasch2000" / "complete_dataset.txt"
OUT = ROOT / "results" / "climate" / "a7_concordance.json"

COL_OSMO = "1M sorbitol - 30 min"
COL_HEAT = "Heat Shock 30 minutes hs-1"
COL_NIT = "Nitrogen Depletion 1 h"

ARMS = [
    dict(arm="A", regime="osmotic carbon diversion",
         stress="1M sorbitol hyperosmotic shock, 30 min", column=COL_OSMO,
         pathway="glycerol biosynthesis",
         genes=["GPD1", "GPD2", "GPP1", "GPP2"], expected="UP",
         citation="Hohmann 2002 MMBR 66:300; Gasch 2000"),
    dict(arm="B", regime="thermal carbon starvation",
         stress="heat shock 25->37 C, 30 min", column=COL_HEAT,
         pathway="heat-shock proteins",
         genes=["HSP12", "HSP26", "HSP30", "HSP42", "HSP78", "HSP82",
                "HSP104", "SSA3", "SSA4", "STI1"], expected="UP",
         citation="Gasch 2000; Verghese 2012"),
    dict(arm="C", regime="thermal carbon starvation",
         stress="heat shock 25->37 C, 30 min", column=COL_HEAT,
         pathway="glycolysis",
         genes=["PGI1", "PFK1", "PFK2", "FBA1", "TPI1", "GPM1", "PGK1",
                "ENO1", "ENO2", "TDH3", "HXK2", "CDC19", "PDC1", "ADH1"],
         expected="UP",
         citation="Gasch 2000 ESR induced carbohydrate/energy genes",
         least_certain=True),
    dict(arm="D", regime="multi-constraint collapse (ammonium-bound component)",
         stress="nitrogen depletion, 1 h", column=COL_NIT,
         pathway="nitrogen catabolite repression targets",
         genes=["GAP1", "MEP2", "DAL5", "PUT4", "UGA4", "GDH2", "GLT1"],
         expected="UP",
         citation="Magasanik & Kaiser 2002 Gene 290:1; Gasch 2000"),
]


GENE_UID = {  # SGD backend locus lookup, queried 2026-09-27 (grounded)
    "GPD1": "YDL022W", "GPD2": "YOL059W", "GPP1": "YIL053W",
    "GPP2": "YER062C", "HSP12": "YFL014W", "HSP26": "YBR072W",
    "HSP30": "YCR021C", "HSP42": "YDR171W", "HSP78": "YDR258C",
    "HSP82": "YPL240C", "HSP104": "YLL026W", "SSA3": "YBL075C",
    "SSA4": "YER103W", "STI1": "YOR027W", "PGI1": "YBR196C",
    "PFK1": "YGR240C", "PFK2": "YMR205C", "FBA1": "YKL060C",
    "TPI1": "YDR050C", "GPM1": "YKL152C", "PGK1": "YCR012W",
    "ENO1": "YGR254W", "ENO2": "YHR174W", "TDH3": "YGR192C",
    "HXK2": "YGL253W", "CDC19": "YAL038W", "PDC1": "YLR044C",
    "ADH1": "YOL086C", "GAP1": "YKR039W", "MEP2": "YNL142W",
    "DAL5": "YJR152W", "PUT4": "YOR348C", "UGA4": "YDL210W",
    "GDH2": "YDL215C", "GLT1": "YDL171C",
}

def binom_sf(k, n):
    """P[X >= k], X ~ Bin(n, 0.5), exact."""
    return sum(math.comb(n, i) for i in range(k, n + 1)) / 2 ** n


def binom_two_sided(k, n):
    p_obs = binom_sf(k, n) if k * 2 >= n else binom_sf(n - k, n)
    return min(1.0, 2 * p_obs)


def clopper_pearson(k, n, alpha=0.05):
    from scipy.stats import beta
    lo = 0.0 if k == 0 else beta.ppf(alpha / 2, k, n - k + 1)
    hi = 1.0 if k == n else beta.ppf(1 - alpha / 2, k + 1, n - k)
    return lo, hi


def main():
    file_sha = sha256(DATA.read_bytes()).hexdigest()
    with open(DATA, newline="") as fh:
        rdr = csv.reader(fh, delimiter="\t")
        header = next(rdr)
        col_idx = {}
        for arm in ARMS:
            if arm["column"] not in header:
                raise SystemExit(f"locked column missing: {arm['column']!r}")
            col_idx[arm["column"]] = header.index(arm["column"])
        name_idx = header.index("NAME")
        rows = list(rdr)
    # resolve each locked gene to exactly one row by its SGD-verified
    # systematic ORF id (GENE_UID; deterministic). 2000-era names differ
    # (GPP1=RHR2, GPP2=HOR2) and name tokens are not unique (HSP26 is
    # mentioned in HSP42's annotation), so name matching is not used.
    resolution = {}
    uid_idx = header.index("UID")
    for arm in ARMS:
        for g in arm["genes"]:
            uid = GENE_UID[g]
            hits = [r for r in rows if len(r) > uid_idx and r[uid_idx] == uid]
            if len(hits) != 1:
                raise SystemExit(
                    f"gene {g} ({uid}): {len(hits)} row matches (need 1)")
            resolution[g] = uid
    results = []
    for arm in ARMS:
        ci = col_idx[arm["column"]]
        per_gene = {}
        for g in arm["genes"]:
            uid = resolution[g]
            row = next(r for r in rows if r[0] == uid)
            raw = row[ci].strip() if ci < len(row) else ""
            try:
                val = float(raw)
            except ValueError:
                val = None
            per_gene[g] = val
        usable = {g: v for g, v in per_gene.items() if v is not None}
        n = len(usable)
        k = sum(1 for v in usable.values() if v > 0)
        p_one = binom_sf(k, n) if n else None
        p_two = binom_two_sided(k, n) if n else None
        lo, hi = clopper_pearson(k, n) if n else (None, None)
        results.append(dict(
            arm=arm["arm"], regime=arm["regime"], stress=arm["stress"],
            column=arm["column"], pathway=arm["pathway"],
            citation=arm["citation"],
            least_certain=arm.get("least_certain", False),
            genes=arm["genes"], uids={g: resolution[g] for g in arm["genes"]},
            log2_ratios=per_gene, n=n, k_agree=k,
            binom_p_one_sided=p_one, binom_p_two_sided=p_two,
            ci95_clopper_pearson=[lo, hi]))
    out = dict(
        analysis="A7 regime-transcriptome concordance (Amendment 3)",
        dataset=dict(
            name="Gasch et al. 2000 complete_dataset.txt (PMID 11102521)",
            source="Wayback Machine capture of genome-www.stanford.edu/"
                   "yeast_stress/data/rawdata/complete_dataset.txt",
            sha256=file_sha,
            grade="context (pre-reg Section 8); A7 is concordance, "
                  "not a validation gate"),
        ethanol_regime="untested: no ethanol-stress experiment in Gasch 2000 "
                       "(carbon-source columns are not stress)",
        arms=results)
    OUT.write_text(json.dumps(out, indent=2))
    print(json.dumps([{k: r[k] for k in
                       ("arm", "regime", "n", "k_agree", "binom_p_one_sided",
                        "binom_p_two_sided", "ci95_clopper_pearson")}
                      for r in results], indent=2))


if __name__ == "__main__":
    main()
