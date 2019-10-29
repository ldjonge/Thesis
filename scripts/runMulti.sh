#!/usr/bin/env bash
old=$(ls multiOutput/raw/ | wc -l)+1
old=$(($old +1))
new=$(($old+100))
for i in $(seq $old $new); do
  python3 -u scripts/newLearning/plotMulti.py;
done
head -1 multiOutput/raw/sim$new.tsv | sed -r "s/^/Run\t/" > multiOutput/summary/summary.tsv
for i in $(seq $old $new); do
  tail -n +2 "multiOutput/raw/sim${i}.tsv" | sed -r "s/^/${i}\t/">> multiOutput/summary/summary.tsv
done
python3 scripts/reshapeSumm.py
