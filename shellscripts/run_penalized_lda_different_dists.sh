#!/bin/sh
# Jake Yeung
# run_penalized_lda_different_dists.sh
# Run penalized LDA for liver-rhythmic genes across distances 
# 2016-04-06

# runscript="/home/yeung/projects/tissue-specificity/scripts/penalized_lda_analysis/multigene_analysis.shellscript.R"
runscript="/home/yeung/projects/tissue-specificity/scripts/penalized_lda_analysis/do_plda_on_posterior_cutoff_0.1.shellscript.R"

[[ ! -e $runscript ]] && echo "$runscript not found, exiting" && exit 1

n=0
maxjobs=6

for dist in 2500 5000 10000; do
	for cutoff in 1.5 2 2.5; do
		echo $dist
		Rscript $runscript $dist $cutoff&
		# limit jobs
		if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
		    wait # wait until all have finished (not optimal, but most times good enough)
		    echo $n wait
		fi
	done
done
wait
