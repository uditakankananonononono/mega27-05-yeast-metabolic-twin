# A6 run log

## 2026-09-27 ~20:30 IST - quantitative overlap complete (honest null)
- Locked script committed pre-run (0a92c26); two pre-outcome bug fixes
  (main-process WT eval, per-regime WT record) committed before results.
- All 4 quantitative comparisons computed on locked universes:
  * ethanol_gse151784: universe 705, model-sensitive 0, dataset hits 604, overlap 0
  * heat_gibney2013: universe 678, model-sensitive 0, dataset hits 536, overlap 0
  * sorbitol_gse59659: universe 1060, model-sensitive 3, dataset hits 141, overlap 0
  * yp20_gse59659: universe 1060, model-sensitive 3, dataset hits 146, overlap 0
- G1 overlap clause: 0/4 significant (BH-FDR 5%) - FAILS, reported plainly
  per Section 10. The model's stress response is dominated by env-imposed
  caps (A4 mechanism), so single-gene deletions produce almost no
  stress-conditional growth sensitivity - consistent with the A5
  all-rel=1.0 negative. The deletion-screen biology (protein synthesis,
  vacuole, mitochondrial genes) lives outside what FBA with env caps can
  see; declared as a model limitation, not hidden.
- Remaining Section-8 items: van Voorst 2008 qualitative gene-set
  concordance (named genes); E-GEOD-1723 and E-MTAB-8245 qualitative
  nutrient checks; held-out narrative checks (no fitting).
