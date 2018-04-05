# Run script: /home/yeung/projects/tissue-specificity/scripts/dhs/dhs_merged_find_peaks.R in batch
# 2015-05-28
# Jake Yeung

beddir="/home/yeung/projects/tissue-specificity/data/beds/merge/ucsc_names"
outmain="/home/yeung/projects/tissue-specificity/data/beds/merge/cutoffs_smoother_gamma_cutoff_high_accuracy_ALL_gamma_inits_scale50_mu3_12_zeros_normalized_scaled_billion"
[[ ! -d $outmain ]] && mkdir $outmain
peakscript="/home/yeung/projects/tissue-specificity/scripts/dhs/dhs_merged_find_peaks.R"

[[ ! -e $peakscript ]] && echo "$peakscript not found, exiting" && exit 1
[[ ! -d $outmain ]] && echo "$outmain not found, exiting" && exit 1
[[ ! -d $beddir ]] && echo "$beddir not found, exiting" && exit 1

[[ ! -d $outmain ]] && mkdir $outmain

for bedin in `ls -d $beddir/*bed`
do
	bedin_full=`readlink -e $bedin`
	bedin_base=$(basename $bedin)
	tissue=${bedin_base%%.*}
	outdir=$outmain/$tissue
	[[ -d $outdir ]] && echo "$outdir found, continuing" && continue
	outbed="$outmain/$tissue/$tissue.dhs.merged.filtered.bed"
	Rscript $peakscript $bedin_full $outdir $outbed&
done
