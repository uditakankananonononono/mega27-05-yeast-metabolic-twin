"""Best third KO on the SDH2+ACH1 design -> results/strain_design_triples.csv"""
import sys, time
sys.path.insert(0, "src")
import pandas as pd
from yeasttwin.media import minimal_glucose
from yeasttwin.model import load_model
from yeasttwin.strain_design import obligatory_production

t0 = time.time()
m0 = minimal_glucose(load_model())
BASE = ("YLL041C", "YBL015W")  # SDH2 + ACH1
b_s, b_g = obligatory_production(m0, BASE, medium=None)
print("SDH2+ACH1 baseline: obligatory=%.4f growth=%.4f" % (b_s, b_g), flush=True)
rows = []
for i, gid in enumerate([g.id for g in m0.genes if g.id not in BASE]):
    s, g = obligatory_production(m0, BASE + (gid,), medium=None)
    rows.append({"gene": gid, "growth": g, "obligatory_succinate": s,
                 "gain_over_pair": s - b_s})
    if (i+1) % 300 == 0:
        print(i+1, "done %.0fs" % (time.time()-t0), flush=True)
df = pd.DataFrame(rows).sort_values("gain_over_pair", ascending=False)
df.to_csv("results/strain_design_triples.csv", index=False)
print(df.head(12).to_string(index=False), flush=True)
print("TOTAL %.0fs" % (time.time()-t0))
