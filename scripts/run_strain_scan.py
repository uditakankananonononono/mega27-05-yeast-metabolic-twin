"""Growth-coupled succinate scan -> results/strain_design_coupled_scan.csv +
best 5th KO on the Raab quadruple background -> results/strain_design_5th_ko.csv"""
import sys, time
sys.path.insert(0, "src")
import pandas as pd
from yeasttwin.media import complete_y7
from yeasttwin.model import load_model
from yeasttwin.strain_design import (find_succinate_exchange,
                                     growth_coupled_scan)

t0 = time.time()
m = load_model()
df = growth_coupled_scan(m)
df.to_csv("results/strain_design_coupled_scan.csv", index=False)
print("coupled scan:", len(df), "genes in %.0fs" % (time.time() - t0), flush=True)

RAAB4 = ("YKL148C", "YLL041C", "YOR136W", "YDL066W")
fifth = growth_coupled_scan(load_model(), background=RAAB4)
fifth.to_csv("results/strain_design_5th_ko.csv", index=False)
# baseline: quadruple alone under the same metric
m2 = complete_y7(load_model())
for bg in RAAB4:
    m2.genes.get_by_id(bg).knock_out()
g4 = m2.slim_optimize(error_value=0.0)
m2.reactions.get_by_id("r_2111").lower_bound = 0.9 * g4
m2.objective = m2.reactions.get_by_id(find_succinate_exchange(m2))
s4 = m2.slim_optimize(error_value=0.0)
print("quadruple baseline: growth=%.4f coupled succinate=%.4f" % (g4, s4), flush=True)
print(fifth.head(10).to_string(index=False), flush=True)
print("TOTAL %.0fs" % (time.time() - t0))
