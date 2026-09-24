#!/bin/bash
cd ~/mega27/item5-yeast
for i in $(seq 1 240); do
  grep -aq EXIT results/cv_run.log && break
  sleep 10
done
PYTHONPATH=src python3 scripts/run_strain_scan.py > results/strain_scan.log 2>&1
echo EXIT_$? >> results/strain_scan.log
