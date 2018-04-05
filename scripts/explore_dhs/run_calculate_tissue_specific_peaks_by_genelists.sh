#!/bin/sh
# Jake Yeung
# run_calculate_tissue_specific_peaks_by_genelists.sh
# Run it 
# 2016-05-09

inscript="/home/yeung/projects/tissue-specificity/scripts/explore_dhs/calculate_liver_specific_peaks_by_genelists.R"
ntrials=10000

n=0
maxjobs=6

# for tiss in Mus Kidney Aorta Hypo BS Cere; do 
for tiss in Liver Lung Mus Kidney Heart Cere; do
	Rscript $inscript $tiss $ntrials&
	if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
		# define maxjobs and n using maxjobsn skeleton
		wait # wait until all have finished (not optimal, but most times good enough)
		echo $n wait
	fi
	ret=$?; [[ $ret -ne 0  ]] && echo "ERROR: script failed" && exit 1
done
wait

