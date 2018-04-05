#!/bin/sh
# Jake Yeung
# copy_motifs_under_peaks_from_vitalit.sh
# to do penalized LDA analysis, we need N matrix of motifs under each peak.
# Because of large file size, we prefilter by gene list (in vitalit example output: 
# /scratch/el/monthly/jyeung/motevo_dhs_outputs/motevo_outputs/liver_kidney_wtko_closestbed_multiple_genes_ncond_modules
# which is generated from 
# /Home/jyeung/projects/tissue_specificity_hogenesch_shellscripts/motevo_dhs_scripts_clean/3-run_motevo/6d-filter_by_gene_list.sh
# 2016-06-26

# Copy output directory into cat'd bedfile to Rstudio

# beddir="/scratch/el/monthly/jyeung/motevo_dhs_outputs/motevo_outputs/liver_kidney_wtko_closestbed_multiple_genes_ncond_modules"  # in vitalit
# outf="/home/yeung/data/tissue_specificity/motevo_dhs/liver_kidney_wtko_modules/all_sites.closest.filter.genelst.dist.50000.long.bed"
beddir=$1
outf=$2
srvr="jyeung@frt.el.vital-it.ch"

ssh $srvr "cat $beddir/*.bed" > $outf
