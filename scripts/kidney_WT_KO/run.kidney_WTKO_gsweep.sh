#!/bin/sh
# Jake Yeung
# run.liver_kidney_cedric.script.sh
# Run nconds
# 2016-06-16

# ncondscript="/home/yeung/projects/tissue-specificity/scripts/liver_kidney_WTKO/liver_kidney_WTKO_nconds.R"
ncondscript="/home/yeung/projects/tissue-specificity/scripts/kidney_WT_KO/nconds_kidney_WT_KO.R"
n=0
maxjobs=9
for g in 101 251 501 751 1001 2501 5001 7501 10001; do
	meth="g=$g"
	Rscript $ncondscript $meth&
	if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
	    # define maxjobs and n using maxjobsn skeleton
	    wait # wait until all have finished (not optimal, but most times good enough)
	    echo $n wait
	fi
done

