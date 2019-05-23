#!/usr/bin/env bash
for i in {1..100}; do
  python3 -u scripts/testPlot.py;
done
scripts/summary.sh
python3 scripts/summPlot.py
