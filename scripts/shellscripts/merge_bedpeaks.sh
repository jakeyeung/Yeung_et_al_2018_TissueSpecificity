# merge_bedpeaks.sh
# 2015-05-28
# Jake Yeung

# After filtering peaks using a cutoff, we merge overlapping or touching peaks together.

filtered_bed_dir="/home/yeung/projects/tissue-specificity/data/beds/merge/cutoffs_manual"
merged_bed_dir="/home/yeung/projects/tissue-specificity/data/beds/merge/cutoffs_manual/merged_peaks"

for f in `ls -d $filtered_bed_dir/*.bed`
do
	bname=$(basename $f)
	tissue=${bname%%.*}
	bedtools merge -c 4 -o sum -i $f > $merged_bed_dir/$tissue.merged_peaks.bed&
done
