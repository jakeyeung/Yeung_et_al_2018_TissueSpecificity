# run_samtools_merge.sh
# Jake Yeung
# 2015-05-29

beddir="/home/yeung/projects/tissue-specificity/data/beds/merge/cutoffs_manual"
outdir="/home/yeung/projects/tissue-specificity/data/beds/merge/cutoffs_manual/merged_peaks"

for f in `ls -d $beddir/*filtered.manual.bed`
do
	bname=$(basename $f)
	tissue=${bname%%.*}
	bedtools merge -c 4 -o sum -i $f > $outdir/"$tissue".merged_peaks.bed&
done

