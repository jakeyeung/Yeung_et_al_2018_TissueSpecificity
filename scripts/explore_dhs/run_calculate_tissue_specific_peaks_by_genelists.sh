#!/bin/sh
# Jake Yeung
# run_calculate_tissue_specific_peaks_by_genelists.sh
# Run it 
# 2016-05-09

inscript="/home/yeung/projects/tissue-specificity/scripts/explore_dhs/calculate_liver_specific_peaks_by_genelists.R"
for tiss in Kidney Aorta Hypo BS Cere; do 
	Rscript calculate_liver_specific_peaks_by_genelists.R $tiss
done

