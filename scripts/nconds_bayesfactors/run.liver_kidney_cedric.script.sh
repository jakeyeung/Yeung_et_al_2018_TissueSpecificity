#!/bin/sh
# Jake Yeung
# run.liver_kidney_cedric.script.sh
# Run nconds
# 2016-06-16

n=0
maxjobs=9
inf="/home/yeung/projects/tissue-specificity/scripts/nconds_bayesfactors/liver_kidney_cedric.script.R"
for g in 10 50 101 251 501 751 1001 5001 10001; do
	meth="g=$g"
	Rscript $inf $meth&
	if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
	    # define maxjobs and n using maxjobsn skeleton
	    wait # wait until all have finished (not optimal, but most times good enough)
	    echo $n wait
	fi

done

