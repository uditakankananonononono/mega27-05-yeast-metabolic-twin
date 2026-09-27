
## Evo stage completion note (2026-09-27 18:5x IST)
- All 3 seeds x 12 gens completed; per-gen screen files present.
- Engineering caveat (no prereg deviation): sandbox timeouts killed the
  runner mid-generation several times. Resume re-screens only unfinished
  designs of the in-flight generation (pop is rebuilt deterministically
  from the previous generation's scores + seed rng, so the bred pop is
  unchanged); completed generations replay from cache. Screen files are
  rewritten deduplicated on replay (duplicate bred designs collapse),
  so per-gen files may show fewer rows than the 24-strong pop even though
  all 24 bred designs were evaluated (duplicates share identical scores).
- Two list-hashability fixes were committed pre-outcome (resume breeding
  normalization, finalist aggregation) - no outcome-affecting change.
- Evo confirm (10 finalists on design150): all rel=1.0, 0 new collapses.
  Evolutionary search over confirmed mods finds nothing above WT - honest
  negative consistent with singles/doubles/scaling screens.

## TEST stage / G4 verdict (2026-09-27 ~20:08 IST)
- 25 confirmed candidates evaluated on test150 (touched once) + tolerance scans.
- Literal locked-rule count: 1/25 passes G4(a)+G4(b)+plausibility
  (YGL248W KO: test rel 1.0000000000000102 vs baseline p95 1.0000000000000018).
- The margin is 8e-15 - inside the declared 1e-14 solver-noise band - and
  all tolerances are identical to WT (37C / 1% EtOH / 0.05 N).
- Honest verdict: G4 FAILS. No candidate improves any objective beyond
  numerical noise. Consistent with singles (1,143), doubles (435),
  scaling (36), evo (3 seeds x 12 gens), all at rel=1.0.
- A4's explanation holds: collapse is driven by env-imposed caps, so
  knockouts/scalings cannot gain. WT tolerance values are mapping
  artifacts by construction (declared, not discoveries).
- filter_b: PENDING - Section-8 datasets not yet acquired (A6).
- Engineering fix pre-outcome: test stage made resume-safe
  (test_progress.jsonl); no locked-rule change.
