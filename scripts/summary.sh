#!/usr/bin/env bash
head -1 output/raw/sim1.tsv | sed -r "s/^/Run\t/" > output/summary/summary.tsv
for i in {101..200}; do
  tail -n +2 "newOutputAgain/raw/sim${i}.tsv" | sed -r "s/^/${i}\t/">> output/summary/summary.tsv
done
