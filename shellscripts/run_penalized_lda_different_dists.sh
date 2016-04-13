#!/bin/sh
# Jake Yeung
# run_penalized_lda_different_dists.sh
# Run penalized LDA for liver-rhythmic genes across distances 
# 2016-04-06

runscript="/home/yeung/projects/tissue-specificity/scripts/penalized_lda_analysis/multigene_analysis.shellscript.R"

[[ ! -e $runscript ]] && echo "$runscript not found, exiting" && exit 1

for dist in 2500 5000; do
	for cutoff in 0.5 2.5 3.5; do
		echo $dist
		Rscript $runscript $dist $cutoff&
	done
done
wait
