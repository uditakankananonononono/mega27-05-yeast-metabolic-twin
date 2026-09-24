"""Full single-KO succinate scan -> results/strain_design_ko_scan.csv"""
import sys, time
sys.path.insert(0, "src")
from yeasttwin.model import load_model
from yeasttwin.strain_design import knockout_scan

t0 = time.time()
df = knockout_scan(load_model())
df.to_csv("results/strain_design_ko_scan.csv", index=False)
print("scanned", len(df), "genes in %.0fs" % (time.time() - t0))
print(df.head(15).to_string(index=False))
