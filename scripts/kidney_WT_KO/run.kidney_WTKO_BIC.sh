#!/bin/sh
# Jake Yeung
# run.liver_kidney_cedric.script.sh
# Run nconds
# 2016-06-16

# ncondscript="/home/yeung/projects/tissue-specificity/scripts/liver_kidney_WTKO/liver_kidney_WTKO_nconds.R"
ncondscript="/home/yeung/projects/tissue-specificity/scripts/kidney_WT_KO/nconds_kidney_WT_KO.R"
n=0
# maxjobs=1
meth="BIC"
Rscript $ncondscript $meth
