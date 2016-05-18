#!/bin/sh
# Jake Yeung
# copy_closestbeds_to_rstudio.sh
# Copy
# 2016-04-29

# indir="/scratch/el/monthly/jyeung/motevo_dhs_outputs/motevo_outputs/closestbed_multiple_genes_Mus_filtered"
tiss="flat"
indir=/scratch/el/monthly/jyeung/motevo_dhs_outputs/motevo_outputs/closestbed_multiple_genes_"$tiss"_filtered
outf="all_sites.closest.filter.$tiss.dist.50000.long.bed"
outdir="/home/yeung/data/tissue_specificity/motevo_dhs/closest_bed_multigene"

ssh jyeung@frt.el.vital-it.ch "cat $indir/*.bed" > $outdir/$outf
