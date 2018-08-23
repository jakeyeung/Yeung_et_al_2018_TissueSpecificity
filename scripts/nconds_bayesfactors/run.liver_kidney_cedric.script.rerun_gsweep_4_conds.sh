#!/bin/sh
# Jake Yeung
# run.liver_kidney_cedric.script.rerun_gsweep_4_conds.sh
# Run it 
# 2018-08-21

runscript="/home/yeung/projects/tissue-specificity/scripts/nconds_bayesfactors/liver_kidney_cedric.script.rerun_gsweep_4_conds.R"

[[ ! -e $runscript ]] && echo "$runscript not found, exiting" && exit 1

for g in 100 250 500 750 1000 1250 1500 2000 2500 5000 10000; do
	Rscript $runscript "g="$g&
done
wait

