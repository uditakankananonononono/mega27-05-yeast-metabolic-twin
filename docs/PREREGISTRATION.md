# Pre-Registration: Virtual Yeast for Climate-Resilient Fermentation

**Program**: MEGA27 item 5 (redirect of mega27-05-yeast-metabolic-twin)
**Repo**: uditakankananonononono/mega27-05-yeast-metabolic-twin (private)
**Locked**: 2026-09-26, Asia/Calcutta. This file is committed and its SHA-256
reported to the program lead BEFORE any outcome data is generated, inspected,
or used to tune any parameter of the pipeline below.
**Owner order**: WhatsApp, 2026-09-26 15:25:51 IST - "delete enviropig and
redirect virtual yeast cell", with the project spec quoted in Section 1.
**Standing constraints**: owner-ordered pause on ChatGPT usage remains in force
(no ChatGPT judge rounds anywhere in this project). Program rules apply: real
open datasets, real code (no stubs/simulations/pseudocode), hermetic pytest
suites, math rigor, honest verdicts, negative results preserved, sandbox
~2 CPU / 1-2 GB RAM.

This document is immutable once committed. Changes are only ever appended as
dated, signed amendments (Section 12). The original text is never edited.

---

## 1. Problem statement (owner spec, verbatim)

"The problem is that yeast used for biofuels, food, and biotechnology can lose
productivity under combinations of heat, ethanol, osmotic stress, and nutrient
limitation. Most studies investigate these stresses separately.

Build a virtual S. cerevisiae cell that integrates metabolism and stress
response. Then computationally expose it to thousands of combined-stress
environments and identify metabolic bottlenecks that cause productivity
collapse.

Then ask: Can we computationally discover metabolic adaptations that allow
yeast to maintain productivity under future environmental conditions?

The impact is concrete: Climate stress -> fermentation failure -> lower
industrial productivity -> identify biological intervention.

Your model could predict specific genes/pathways whose modification should
make yeast more resilient. You could validate those predictions using existing
experimental datasets rather than performing genetic engineering yourself.

Even better: target a specific industry. For example, bioethanol production.
Instead of maximizing growth, optimize: ethanol yield + temperature tolerance
+ ethanol tolerance + low nutrient requirements.

Your virtual cell could search for combinations of metabolic changes that
satisfy all four simultaneously. That gives you a potential practical outcome:
A computationally predicted yeast engineering strategy for maintaining
bioethanol production under heat and ethanol stress."

Framing rule (owner spec): the project is the biological/industrial problem.
The virtual cell is the tool, not the deliverable. ISEF framing follows the
problem, not the model.

## 2. Redirect relationship to the closed item-5 direction

The closed direction (essentiality prediction: GCN/CNN vs FBA benchmark;
isozyme over-rescue discovery; sdh2/ach1 succinate strain design) is finished
and archived under legacy/. Its gates are closed and none of its outcome
credit carries over. What is reused is infrastructure only: the in-repo
yeast-GEM SBML and COBRApy loading/media/FBA harness, the hermetic test
scaffold, the dataset-manifest/provenance discipline, the Times/lualatex
paper pipeline, and the shipped-CLI pattern. Every scientific claim in this
redirect is new, locked here, and evaluated only against data acquired and
analyzed after this lock.

## 3. Research question and hypotheses

**Question**: Can we computationally discover metabolic adaptations that allow
S. cerevisiae to maintain bioethanol productivity under combined heat,
ethanol, osmotic, and nutrient-limitation stress - and validate those
predictions against existing experimental datasets?

- **H1 (non-additivity)**: combined-stress productivity collapse is not
  predictable from single-stress dose-responses; the deviation is systematic,
  not noise.
- **H2 (bottlenecks)**: productivity collapse under combined stress is
  driven by a small, identifiable set of metabolic bottlenecks (cofactor
  balance, osmoregulatory carbon drain, membrane/lipid demand, translation-
  capacity precursors), discoverable from shadow-price/FVA structure.
- **H3 (engineering)**: there exist small (<=3 modification) in-silico strain
  designs that Pareto-dominate wild type across all four locked objectives
  (Section 6) in held-out environments, with at least one implicated gene
  supported by concordant evidence in at least one locked experimental
  dataset (Section 8).

Null outcomes for any hypothesis are reportable results. Baselines and
non-beats are reported plainly (Section 10).

## 4. Model

- **Base**: the in-repo consensus yeast-GEM SBML (4,105 reactions / 2,748
  metabolites / 1,143 genes; upstream: SysBioChalmers yeast-GEM, releases
  8.7.x, https://github.com/SysBioChalmers/yeast-GEM). The exact model file
  is pinned by SHA-256 in the dataset manifest at first commit of this phase.
- **Stress-response integration layer** (the new build): quantitative,
  source-cited mappings from each stressor axis to model constraints:
  - Heat: temperature-dependent enzyme turnover/deactivation scaling and
    increased non-growth-associated maintenance (NGAM), parameterized from
    published growth-rate-vs-temperature curves for S. cerevisiae; declared
    Arrhenius-style functional form with locked parameters before runs.
  - Ethanol: published inhibition kinetics on growth and fermentation
    (e.g., ATCC 4126 batch study, Biotechnol. Bioeng. 27:3, doi
    10.1002/bit.260270311; biomass-yield decline 0.156 -> 0.026 over 0-107
    g/L ethanol, doi 10.1002/bit.260400213), mapped to uptake/maintenance
    and membrane-transport constraints.
  - Osmotic: water-activity/osmolarity mapped to obligatory glycerol-
    production demand (osmoregulatory carbon drain) and growth-rate
    penalty, parameterized from published osmotic-stress physiology.
  - Nutrient limitation: uptake bounds for N, P, S and vitamin pools, with
    limitation levels defined as locked fractions of the standard medium.
- **Enzyme constraints**: full GECKO-style enzyme-constrained model only if
  it runs within sandbox memory; otherwise the declared fallback is
  selective enzyme constraints on stress-relevant subsystems (glycolysis,
  TCA, glycerol, lipid/membrane, ergosterol, trehalose, amino-acid pools),
  with the selection rule locked in code before runs. The fallback is a
  declared limitation, not hidden.
- Every stress-to-constraint parameter carries a citation to a public
  dataset or published value. Where no published value exists, the
  parameter is marked ASSUMED with a pre-registered sensitivity sweep
  (+-50%), and is never fitted to any validation outcome.

## 5. Environment generator

Four axes, locked ranges justified from industrial fermentation literature
(very-high-gravity bioethanol review, Fermentation 7(1):38; multiple-stress
tolerance panel, AMB Expr. 6:87 doi 10.1186/s13568-016-0285-x):

| Axis | Range | Resolution |
|---|---|---|
| Temperature | 28-42 C | 1 C steps (15 levels) |
| Ethanol (v/v) | 0-14% | 1% steps (15 levels) |
| Osmotic (sorbitol-equivalent) | 0-1.5 M | 0.1 M steps (16 levels) |
| Nitrogen availability | 5-100% of standard | 5% steps (20 levels) |

Full grid = 72,000 environments. Sampling: full grid for single- and
pairwise-stress maps; Latin-hypercube subsample of 5,000 four-axis
combinations (seed 20260926) for the combined-stress collapse surface. Grid
split into DESIGN and TEST regions before any optimization: stratified 70/30
split, seed 424242, committed in code. All optimization and tuning happens
on DESIGN only; TEST is touched once, for final claims.

## 6. Objectives (locked definitions)

Per environment, on the constrained model:

1. **Ethanol yield**: max ethanol secretion flux subject to environment
   constraints and locked minimal viability (ATP maintenance feasible),
   reported as mol EtOH / mol glucose and g/g.
2. **Temperature tolerance**: highest grid temperature at which ethanol
   yield stays >= 80% of the unstressed reference yield (reference: 30 C,
   0% EtOH, 0 M osmotic, 100% N).
3. **Ethanol tolerance**: highest grid ethanol concentration at which the
   same 80% threshold holds.
4. **Nutrient requirement**: lowest grid nitrogen level at which the same
   80% threshold holds.

Biomass maximization is NOT an objective (owner spec). Viability is a
constraint, not a target. An environment is "collapsed" for a genotype when
max ethanol flux is below 20% of the reference value or infeasible.

## 7. Analyses (locked)

- **A1**: single-stress dose-responses for all four axes (WT).
- **A2**: combined-stress collapse surface on the 5,000-environment LHS
  sample (WT).
- **A3 (H1)**: non-additivity test. Independence prediction per environment
  = product of single-stress relative yields (Bliss-style). Deviation =
  simulated combined relative yield minus independence prediction. Locked
  metric: fraction of collapsed-region environments with deviation < 0 and
  bootstrap 95% CI over declared parameter uncertainty excluding 0.
- **A4 (H2)**: bottleneck identification. Shadow prices on uptake/maintenance
  constraints at the collapse boundary, FVA on central carbon and cofactor
  reactions, and conditional essentiality per environment. Bottlenecks
  clustered across the collapse region; cluster membership counts reported.
- **A5 (H3)**: strain-design search. Modification space: single and double
  gene deletions (model genes), reaction-bound scaling as an
  up/down-regulation proxy (locked fold-changes: 0.25x, 0.5x, 2x, 4x), and
  GPR-consistent cofactor swaps. Budget: <=3 modifications per design.
  Search: seeded evolutionary multi-objective (Pareto front over the four
  Section-6 objectives) on the DESIGN region, seeds 42, 123, 2026. Final
  claims only on the TEST region.
- **A6 (validation)**: comparisons against the locked datasets of Section 8,
  using only the locked metrics of Section 8.

## 8. Locked validation datasets and metrics (existing experiments only; no wet lab)

All accession-level; acquired with hashes recorded in the dataset manifest
(same discipline as the closed direction). Model parameters are never fitted
to these outcomes. Any dataset added later is EXPLORATORY and labeled so.

| Stressor | Dataset | Accession / source | Locked metric |
|---|---|---|---|
| Ethanol | Deletion-library fitness under acute lethal ethanol | GEO GSE151784 (SGD) | Stress-conditional sensitivity overlap: Fisher OR vs random, BH-FDR 5%; MCC with 95% CI |
| Ethanol | Deletion-set ethanol tolerance screen (limited aeration) | van Voorst et al., Am. J. Enol. Vitic. 59(4):401 | Qualitative gene-set concordance (named genes only) |
| Heat | Bar-seq deletion-pool heat-shock survival screen | Gibney et al., PNAS 2013, doi 10.1073/pnas.1318100110 | Same overlap metric as GSE151784 |
| Osmotic | Deletion fitness under 20% glucose (YP20) and 1M sorbitol | Yoshikawa et al. 2009, GEO GSE59659 (SGD S000128535) | Same overlap metric |
| Nutrient | C/N/P/S-limited chemostat reference physiology | ArrayExpress E-GEOD-1723 | Qualitative direction check on uptake/shadow-price signs |
| Nutrient | N-limitation metabolic/translational capacity reserves | ArrayExpress E-MTAB-8245 | Qualitative: are predicted spare-capacity bottlenecks concordant |
| Context only | Environmental stress response expression program | Gasch et al. 2000 | Context, never validation |

Held-out qualitative check: published industrial/evolved multiple-stress-
tolerant strain reports (AMB Expr. 6:87 panel; VHG literature) are used only
as narrative held-out checks on final designs, never fitted to.

## 9. Success gates (binary, reported plainly)

- **G1 (model validity)**: WT dose-responses qualitatively reproduce known
  stress physiology (monotonic productivity decline per axis; heat and
  ethanol interaction region consistent with industrial literature) AND
  predicted stress-conditional sensitivity significantly overlaps screen
  hits (Fisher OR > 1, BH-FDR 5%) in at least 2 of the 4 quantitative
  stressor datasets of Section 8.
- **G2 (H1)**: non-additivity confirmed in >= 25% of collapsed-region
  environments with bootstrap CI excluding 0; otherwise the null is
  reported.
- **G3 (H2)**: at least one bottleneck cluster supported by both
  shadow-price and FVA evidence across >= 10% of collapsed environments.
- **G4 (H3)**: at least one <=3-modification design that (a) Pareto-
  dominates WT on all four objectives in the TEST region, (b) beats >= 95%
  of matched-budget random designs (>= 30 seeds), and (c) has at least one
  implicated gene with concordant evidence in at least one locked Section-8
  dataset.
- Any gate failed: reported as failed, with full metrics. Negative results
  are preserved in the paper, then the program's pivot rule applies.

## 10. Honest-negative register

Baselines for every claim: WT, single-stress optima, matched-budget random
designs. Non-beats are reported in the same tables as beats, with the same
statistics. No gate, claim, or table may report a beat without its baseline
in the same artifact.

## 11. Declared assumptions and risks

- FBA is steady-state; no dynamics, no gene-regulatory kinetics. This
  bounds all claims to metabolic feasibility.
- Stress-to-constraint parameters derive from different strain backgrounds
  (S288C-class vs industrial); declared per parameter, sensitivity swept.
- Enzyme-constraint fallback (Section 4) limits proteome-allocation claims.
- The independence model in A3 is a locked null, not a biological claim.
- Literature ethanol-inhibition kinetics are batch-derived; mapping to
  constraint-based form is an approximation, declared in the paper.

## 12. Amendment policy

Amendments are appended below with date, reason, and SHA-256 of the amended
file. The original text above is never modified. Outcome data generated
before an amendment is never re-analyzed under the amendment and presented
as pre-registered; it is labeled exploratory.

## 13. Compliance

All datasets are public (GEO, ArrayExpress, SGD, journal supplements);
usage terms checked at acquisition and recorded in the dataset manifest.
No ChatGPT usage anywhere in this project (owner order). No wet-lab claims:
every prediction is computational, validated only against existing
experimental data.

---

## Amendment 1 - 2026-09-26 (appended before the first outcome run)

Precise implementation rules for the Section-4 stress mappings and the
Section-6 viability constraint. Locked before any outcome data was
generated; the original text above is unchanged.

1. **Heat mapping**: f_T (CTMI) scales glycolytic capacity through the
   glucose uptake lower bound (base -20 mmol/gDW/h in the complete_y7
   medium) and, together with f_E, caps growth at f_T x f_E x reference
   max. Rationale: temperature-dependent enzyme turnover/Vmax decline.
2. **Ethanol mapping**: f_E (Levenspiel-type) caps the ethanol exchange
   upper bound at f_E x reference max (product inhibition of the
   fermentation rate, as measured by the cited kinetics literature) and,
   with f_T, caps growth.
3. **Maintenance scaling**: NGAM(env) = 0.7 / max(f_T x f_E, 0.01)
   mmol ATP/gDW/h (stress raises maintenance energy; ASSUMED functional
   form, swept +-50% per Section 4). Declared empirical note: under the
   locked microaerobic base, NGAM rises alone do not reduce max ethanol
   (respiration of non-sugar carbon absorbs it); the productive couplings
   are rules 1 and 2. This was discovered in pre-run implementation
   testing, before any A1-A6 outcome run.
4. **Viability floor**: an environment is nonviable when
   f_T x f_E < 0.01 (growth/no-growth boundary convention from predictive
   microbiology). Nonviable environments are collapsed by definition. This
   refines the Section-6 "ATP maintenance feasible" clause; it does not
   replace it.
5. **Base-medium correction**: standard ammonium availability in the
   complete_y7 base is the benchmark medium's free ammonium exchange
   (lb -1000); nitrogen fractions of Section 5 scale that bound.
6. **Base oxygen**: microaerobic lb -0.5. Fully anaerobic growth is
   infeasible in this medium (verified pre-run); industrial bioethanol
   fermentation is microaerobic. Documented modeling choice.

SHA-256 of this amended file is recorded in the commit message and
reported to the program lead.

---

## Amendment 2 - 2026-09-26 (appended before the sweep and A4 runs)

1. **Sensitivity-sweep design** (Section 4 rule made precise): one-at-a-time
   sweeps of the three ASSUMED parameters - k_gly {0.5, 2.0} (nominal 1.0),
   ethanol n {0.75, 2.25} (nominal 1.5), T_max {40, 44} (nominal 42) - plus
   the nominal set: 7 parameter sets total. Each set is evaluated on the
   A1 grid curves and a seeded 1,000-environment subsample of the locked
   5,000-environment LHS sample (subsample = every 5th environment by
   sample index, deterministic). The full-5,000 nominal-parameter A3 result
   already reported stands as the nominal verdict; sweeps quantify
   parameter uncertainty around it and are never used to retune nominals.
2. **A4 bounded protocol**: binding-constraint identification at every
   collapsed environment of the subsample (which of: glucose cap, ethanol
   cap, growth cap, NGAM, glycerol drain, ammonium bound binds at the
   optimum, with dual values where the solver provides them), clustered
   across the collapsed set. FVA is run on a declared subset: central
   carbon + fermentation + glycerol + ammonium reactions (locked list in
   code), for 50 seeded collapsed environments (seed 777). G3 unchanged:
   at least one bottleneck cluster supported by both binding/dual evidence
   and FVA across >= 10% of collapsed environments.

## Amendment 3 - 2026-09-26 (appended after judge round 1, before any A5/A7/A8 outcome data)

Prompted by judge round 1 (docs/JUDGE_ROUNDS.md) and folded back under the
novelty rule. Locked BEFORE the analyses it names are run.

1. **Claim reframing.** The central claim is downgraded from "combined
   stress creates universal nonlinear productivity collapse" (G2 framing,
   which FAILED and is not parameter-robust) to: "combined climate stresses
   drive yeast into discrete, condition-dependent metabolic failure regimes
   characterized by competing resource-allocation demands." G2 remains
   reported as failed; no re-gating of G2 is attempted.
2. **A7 (new, runs BEFORE A5): stress-regime topology + conservation.**
   Regimes are defined from the A4 binding-constraint vectors over the
   locked 5,000-environment surface (regime = dominant binding-constraint
   cluster, v1 labels locked in results/climate/a4_regimes.json).
   Validation: concordance between model regimes and published single-stress
   expression programs (glycerol GPD1/GPD2/GPP1/GPP2, heat-shock HSP and
   glycolytic genes, NCR nitrogen genes) from public transcriptome datasets
   (Gasch 2000 is context; the locked validation datasets of Section 8 stay
   as they are). Concordance metric locked: direction-of-change agreement of
   the regime's signature pathway against the matching stress transcriptome,
   reported as a binomial test vs 0.5 with 95% CI. No parameter is refit.
3. **A8 (new): minimal dynamic arm.** A static-optimization dynamic-FBA
   heat-wave trajectory (30 C baseline, 40 C ramp, 42 C peak, 30 C recovery)
   on the locked medium, tracking ethanol, glycerol, growth, and maintenance
   burden over time. Claim limited to: cumulative-exposure threshold
   behavior vs the static A2 surface. No new parameters beyond the locked
   stress layers; trajectory grid locked in code before the run.
4. **A5 reorder + plausibility filter.** A5 strain design now runs after A7
   and adds a locked biological-plausibility screen before any design is
   reported: (a) not essential under the target regime per the model's own
   essentiality call; (b) deletion-screen concordance where the Section 8
   datasets cover the gene; (c) reported with the regime it targets. Designs
   failing the screen are reported as screened-out, never as candidates.
