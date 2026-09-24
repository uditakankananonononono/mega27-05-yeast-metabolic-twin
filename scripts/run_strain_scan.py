"""Full strain-design scan -> results/strain_design_ko_scan.csv +
5th-KO-on-quadruple-background -> results/strain_design_5th_ko.csv"""
import sys, time
sys.path.insert(0, "src")
import pandas as pd
from yeasttwin.media import complete_y7
from yeasttwin.model import load_model
from yeasttwin.strain_design import find_succinate_exchange, knockout_scan

t0 = time.time()
m = load_model()
df = knockout_scan(m)
df.to_csv("results/strain_design_ko_scan.csv", index=False)
print("scanned", len(df), "genes in %.0fs" % (time.time() - t0), flush=True)

# validated background: Raab 2010 quadruple KO
RAAB4 = {"YKL148C": "SDH1", "YLL041C": "SDH2", "YOR136W": "IDH1", "YDL066W": "IDP1"}
m2 = complete_y7(load_model())
wt = m2.slim_optimize()
for orf in RAAB4:
    m2.genes.get_by_id(orf).knock_out()
quad_growth = m2.slim_optimize(error_value=0.0)
ex = find_succinate_exchange(m2)
m2.reactions.get_by_id("r_2111").lower_bound = 0.1 * wt
m2.objective = m2.reactions.get_by_id(ex)
quad_succ = m2.slim_optimize(error_value=0.0)
print("quadruple background: growth=%.4f max succinate=%.4f" % (quad_growth, quad_succ), flush=True)

rows = []
for gid in [g.id for g in m2.genes if g.id not in RAAB4]:
    with m2:
        try:
            m2.genes.get_by_id(gid).knock_out()
        except KeyError:
            continue
        g = m2.slim_optimize(error_value=0.0)
        if g < 0.1 * wt:
            continue
        m2.reactions.get_by_id("r_2111").lower_bound = 0.1 * wt
        m2.objective = m2.reactions.get_by_id(ex)
        s = m2.slim_optimize(error_value=0.0)
        rows.append({"gene": gid, "growth": g, "succinate_flux": max(s, 0.0)})
fifth = pd.DataFrame(rows).sort_values("succinate_flux", ascending=False)
fifth["baseline_quadruple_succinate"] = quad_succ
fifth["gain_over_quadruple"] = fifth["succinate_flux"] - quad_succ
fifth.to_csv("results/strain_design_5th_ko.csv", index=False)
print("5th-KO scan done:", fifth.head(10).to_string(index=False), flush=True)
print("TOTAL %.0fs" % (time.time() - t0))
