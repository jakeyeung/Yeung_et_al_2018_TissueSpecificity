#!/bin/sh
# Jake Yeung
# copy_to_vitalit.sh
# Copy gene lists to vitalit 
# 2016-06-26

remotehost=jyeung@frt.el.vital-it.ch
# indir="/home/yeung/projects/tissue-specificity/data/gene_lists/liver_kidney_wtko_modules"
# outf="/scratch/el/monthly/jyeung/motevo_dhs_outputs/gene_lists/liver_kidney_wtko/liver_kidney_wtmodules_bic_g1001.txt"
# indir="/home/yeung/projects/tissue-specificity/data/gene_lists/liver_kidney_wtko_modules/by_g1001"
# outf="/scratch/el/monthly/jyeung/motevo_dhs_outputs/gene_lists/liver_kidney_wtko/liver_kidney_wtmodules_g1001.txt"
# indir="/home/yeung/projects/tissue-specificity/data/gene_lists/liver_kidney_wtko_modules/by_bic"
# outf="/scratch/el/monthly/jyeung/motevo_dhs_outputs/gene_lists/liver_kidney_wtko/liver_kidney_wtmodules_bic.txt"
indir=$1
outf=$2

[[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1

# cat $indir/*.list 
# cat $indir/*.list | ssh $remotehost "cat - > /scratch/el/monthly/jyeung/motevo_dhs_outputs/gene_lists/liver_kidney_wtko/liver_kidney_wtmodules_bic_g1001.txt"
cat $indir/*.list | ssh $remotehost "cat - | sort | uniq > $outf"
