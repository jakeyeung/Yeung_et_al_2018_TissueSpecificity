#!/bin/sh
# Jake Yeung
# run_nconds.sh
# Run nconds 11 tissues shared parameters
# 2015-10-08

datlongdir="/home/yeung/projects/tissue-specificity/data/nconds2/datlong_11_tiss_by_gene"
chunkdir="/home/yeung/projects/tissue-specificity/data/nconds2/mats_11_tiss_max_3_chunks_20000"
outmain="/home/yeung/projects/tissue-specificity/data/nconds2/fits_11_tiss_chunks.max_3_test.20000"
runscript="/home/yeung/projects/tissue-specificity/scripts/nconds/load_chunk_datgene_run_nconds.R"

[[ ! -e $runscript ]] && echo "$runscript not found, exiting" && exit 1
[[ ! -d $chunkdir ]] && echo "$chunkdir not found, exiting" && exit 1
[[ -d $outmain ]] && echo "$outmain" must be empty folder && exit 1
[[ ! -d $outmain ]] && mkdir $outmain

n=0
maxjobs=15

for genepath in `ls -d $datlongdir/*.Robj`; do
	Rscript $runscript $genepath $chunkdir $outmain&
	# limit jobs
	if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
		wait # wait until all have finished (not optimal, but most times good enough)
		echo $n wait
	fi	
done
wait
