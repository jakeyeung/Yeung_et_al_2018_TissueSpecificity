#!/bin/sh
# Jake Yeung
# run.goterms.all.sh
#  
# 2016-11-06

indir="/home/yeung/projects/tissue-specificity/scripts/gene_enrichment_analysis"
for module in kidneyWT liverKO liverWT; do
	fname=`echo "go_terms_analysis_different_phases."$module".all.R"`
	inf=$indir/$fname
	[[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1
	Rscript $inf
done       
