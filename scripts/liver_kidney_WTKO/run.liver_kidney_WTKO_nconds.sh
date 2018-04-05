#!/bin/sh
# Jake Yeung
# run.liver_kidney_WTKO_nconds.sh
#  
# 2016-06-23

ncondscript="/home/yeung/projects/tissue-specificity/scripts/liver_kidney_WTKO/liver_kidney_WTKO_nconds.R"

n=0
maxjobs=2
for meth in BIC zf; do
	Rscript $ncondscript $meth&
	if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
	    # define maxjobs and n using maxjobsn skeleton
	    wait # wait until all have finished (not optimal, but most times good enough)
	    echo $n wait
	fi
done
