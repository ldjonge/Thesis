#!/usr/bin/env bash
head -1 output/raw/sim1.tsv | sed -r "s/^/Run\t/" > summary.tsv
for i in {1..100}; do
  tail -n +2 "output/raw/sim${i}.tsv" | sed -r "s/^/${i}\t/">> summary.tsv
done
