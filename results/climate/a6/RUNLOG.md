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

## van Voorst/AJEV qualitative concordance (2026-09-27 ~20:35 IST)
- Prereg attribution note: Am. J. Enol. Vitic. 59(4):401 is Renuka Kumar
  et al. (Bisson lab) 2008, not van Voorst; the locked accession was used.
- Named genes from the paper (abstract; full text paywalled): CLC1, GSH1,
  UME6, TDP3, VPS24 (persistent weak growth under 5% ethanol), ADH1.
- SGD-verified systematic IDs: CLC1=YGR167W, GSH1=YJL101C, UME6=YDR207C,
  VPS24=YKL041W, ADH1=YOL086C; TDP3 unresolvable via SGD API (recorded).
- Concordance (model ethanol exemplar 8%):
  * ADH1: in model, KO grows = WT under both ref and ethanol (sens=1.0;
    isozyme redundancy ADH2/3/4/5) -> model NOT sensitive, paper sensitive:
    discordant.
  * GSH1: KO unconditionally zero-growth in model (excluded per locked
    rule); likely a model artifact (no GSH salvage/uptake); declared.
  * CLC1, UME6, VPS24: not in the metabolic model (channel/TF/ESCRT).
  * TDP3: unresolved.
- Result: 0 of 1 evaluable named genes concordant - honest null, same
  declared limitation as the quantitative overlap.

## Nutrient qualitative checks - LOCKED DESIGN (pre-run)
- E-GEOD-1723 direction check (locked rule): under locked N-, P-, S-
  limited regimes (25% of standard), the model's limiting-nutrient uptake
  must saturate its bound (uptake sign/shadow price positive on the
  limiting nutrient) and non-limiting uptakes must not; compare sign
  pattern to published C/N/P/S chemostat physiology qualitatively.
- E-MTAB-8245 spare-capacity check (locked rule): identify predicted
  bottleneck reactions under N-limitation via shadow prices on N uptake
  and glycolytic/TCA flux distributions; qualitatively compare against
  the published finding of translational/metabolic capacity reserves
  under N limitation (documented as narrative concordance, no fitting).

## Nutrient checks result (2026-09-27 ~20:40 IST) - MATERIAL LIMITATION FOUND
- E-GEOD-1723 direction check: NOT EVALUABLE. Under locked N-limitation
  (nitrogen_frac=0.25) all fluxes and growth are bit-identical to ref:
  the ammonium lever (r_1654 lb = -1000*frac) is INERT because the locked
  rich base medium supplies nitrogen via amino-acid uptakes (glutamate,
  glutamine, etc., several at their -0.5 bounds); ammonium is EXCRETED
  (flux +3.76), never consumed, so restricting its uptake binds nothing.
- Consequence (declared, not silently fixed): the entire nitrogen axis of
  the locked environment grid is effectively unperturbed in A2/A4/A5.
  The previously declared "N-tol saturates at 5% grid edge" mapping
  artifact is hereby root-caused: the N lever does nothing in this
  medium. Retrofitting the medium/lever post-outcome would violate the
  locked-mapping discipline; flagged to parent as material for the
  courier paste and the redirect paper (limitation + possible future
  amendment: scale amino-acid N or use defined minimal medium).
- E-MTAB-8245 spare-capacity check: also not evaluable for the same
  reason (no N-limited state exists in the model).
- Script note: glucose exchange corrected to r_1714, sulphate r_2060.
