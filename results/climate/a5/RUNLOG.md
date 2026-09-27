
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
