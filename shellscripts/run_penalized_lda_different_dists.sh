#!/bin/sh
# Jake Yeung
# run_penalized_lda_different_dists.sh
# Run penalized LDA for liver-rhythmic genes across distances 
# 2016-04-06

runscript="/home/yeung/projects/tissue-specificity/scripts/penalized_lda_analysis/multigene_analysis.shellscript.R"

[[ ! -e $runscript ]] && echo "$runscript not found, exiting" && exit 1

for dist in 1501 2501 5001 10001 25001; do
	echo $dist
	Rscript $runscript $dist
done
