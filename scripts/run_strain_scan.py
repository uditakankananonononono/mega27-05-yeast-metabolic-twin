"""Obligatory-production (growth-coupling) scans on minimal glucose medium:
1. all single KOs -> results/strain_design_obligatory_singles.csv
2. 5th KO on the Raab quadruple background -> results/strain_design_5th_ko.csv
3. all pairs among a TCA/glyoxylate target set -> results/strain_design_pairs.csv
"""
import itertools, sys, time
sys.path.insert(0, "src")
import pandas as pd
from yeasttwin.model import load_model
from yeasttwin.strain_design import obligatory_production

t0 = time.time()
RAAB4 = ("YKL148C", "YLL041C", "YNL037C", "YDL066W")  # SDH1 SDH2 IDH1 IDP1 (verified ORFs)

# 1. singles (checkpointed); reuse ONE minimal-medium model, no copies
from yeasttwin.media import minimal_glucose
m0 = minimal_glucose(load_model())
MED = None
rows = []
for i, gid in enumerate([g.id for g in m0.genes]):
    s, g = obligatory_production(m0, (gid,), medium=MED)
    rows.append({"gene": gid, "growth": g, "obligatory_succinate": s})
    if (i + 1) % 100 == 0:
        pd.DataFrame(rows).to_csv("results/strain_design_obligatory_singles.csv", index=False)
        print("singles %d/1143 %.0fs" % (i+1, time.time()-t0), flush=True)
df = pd.DataFrame(rows).sort_values("obligatory_succinate", ascending=False)
df.to_csv("results/strain_design_obligatory_singles.csv", index=False)
print("singles done %.0fs; >0: %d" % (time.time()-t0, (df.obligatory_succinate>1e-6).sum()), flush=True)

# 2. 5th KO on quadruple background
q_s, q_g = obligatory_production(m0, RAAB4, medium=MED)
print("quadruple baseline: obligatory=%.4f growth=%.4f" % (q_s, q_g), flush=True)
rows = []
for i, gid in enumerate([g.id for g in m0.genes if g.id not in RAAB4]):
    s, g = obligatory_production(m0, RAAB4 + (gid,), medium=MED)
    rows.append({"gene": gid, "growth": g, "obligatory_succinate": s,
                 "gain_over_quadruple": s - q_s})
    if (i + 1) % 200 == 0:
        pd.DataFrame(rows).to_csv("results/strain_design_5th_ko.csv", index=False)
        print("5th %d/1139 %.0fs" % (i+1, time.time()-t0), flush=True)
fifth = pd.DataFrame(rows).sort_values("gain_over_quadruple", ascending=False)
fifth.to_csv("results/strain_design_5th_ko.csv", index=False)
print("5th-KO done %.0fs" % (time.time()-t0), flush=True)
print(fifth.head(10).to_string(index=False), flush=True)

# 3. pairs among TCA/glyoxylate target set (named by ORF)
TARGETS = ["YKL148C","YLL041C","YNL037C","YDL066W","YOR136W",  # SDH1 SDH2 IDH1 IDP1 IDH2
           "YLR174W","YER065C","YNL117W","YCR005C","YOL126C",   # IDP2 ICL1 MLS1 CIT2 MDH2
           "YDL078C","YPL262W","YOR142W","YGR244C","YIL125W",   # MDH3 FUM1 LSC1 LSC2 KGD1
           "YDR148C","YBL015W","YKR097W","YKL029C","YNR001C",   # KGD2 ACH1 PCK1 MAE1 CIT1
           "YLR304C","YKL120W","YJR095W"]                       # ACO1 OAC1 SFC1 (ORFs verified vs swissprot aliases)
rows = []
for a, b in itertools.combinations(TARGETS, 2):
    s, g = obligatory_production(m0, (a, b), medium=MED)
    rows.append({"ko1": a, "ko2": b, "growth": g, "obligatory_succinate": s})
pairs = pd.DataFrame(rows).sort_values("obligatory_succinate", ascending=False)
pairs.to_csv("results/strain_design_pairs.csv", index=False)
print("pairs done %.0fs" % (time.time()-t0), flush=True)
print(pairs.head(10).to_string(index=False), flush=True)
print("TOTAL %.0fs" % (time.time()-t0))
