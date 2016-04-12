# Jake Yeung
# 19-03-2015
# batch_run_elastic_net.genomewide.sh

sitecounts_dir=/home/yeung/projects/tissue-specificity/site_count_matrices/GLM
output_dir=/home/yeung/projects/tissue-specificity/plots/elastic_net/genome_wide/GLM
rscript=/home/yeung/projects/tissue-specificity/scripts/run_elastic_net.genomewide.R

for d in Sobel.33 SwissRegulon
	do
	for subd in `ls  $sitecounts_dir/$d/`;
		do
		infile=$sitecounts_dir/$d/$subd/sitecount_matrix
		outfile=$output_dir/$d/$subd
		# mkdir $outfile
		for a in 0.005 0.05 0.1 0.5 1
			do
			Rscript $rscript $infile $outfile $a
		done 
	done
done
