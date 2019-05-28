#!/usr/bin/env bash
for i in {101..200}; do
  python3 -u scripts/plotMulti.py;
done
head -1 multiOutput/raw/sim101.tsv | sed -r "s/^/Run\t/" > multiOutput/summary/summary.tsv
for i in {101..200}; do
  tail -n +2 "multiOutput/raw/sim${i}.tsv" | sed -r "s/^/${i}\t/">> multiOutput/summary/summary.tsv
done
python3 scripts/summPlot.py multi
