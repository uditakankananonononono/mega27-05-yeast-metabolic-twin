"""Envelope-shift succinate scans -> results/strain_design_shift_scan.csv
and 5th-KO on Raab quadruple background -> results/strain_design_5th_ko.csv"""
import sys, time
sys.path.insert(0, "src")
from yeasttwin.model import load_model
from yeasttwin.strain_design import envelope_shift_scan

t0 = time.time()
df, wt_val = envelope_shift_scan(load_model())
df.to_csv("results/strain_design_shift_scan.csv", index=False)
print("WT succinate at 50%% growth floor: %.4f" % wt_val, flush=True)
print("shift scan:", len(df), "genes in %.0fs" % (time.time() - t0), flush=True)

RAAB4 = ("YKL148C", "YLL041C", "YOR136W", "YDL066W")
fifth, _ = envelope_shift_scan(load_model(), background=RAAB4)
fifth.to_csv("results/strain_design_5th_ko.csv", index=False)
print("5th-KO scan done", flush=True)
print("TOTAL %.0fs" % (time.time() - t0))
