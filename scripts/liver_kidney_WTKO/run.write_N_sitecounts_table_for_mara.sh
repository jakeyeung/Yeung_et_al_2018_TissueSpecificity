#!/bin/sh
# Jake Yeung
# run.write_N_sitecounts_table_for_mara.sh
# Run it with parameters 
# 2016-07-10

writescript="/home/yeung/projects/tissue-specificity/scripts/liver_kidney_WTKO/write_N_sitecounts_table_for_mara.R"
[[ ! -e $writescript ]] && echo "$writescript not found, exiting" && exit 1

n=0
maxjobs=6

# dist=40000
# cutoff=3
model="Liver_SV129"
for method in "g=1001" "BIC"; do
	for dist in 10000 20000 40000; do
		for cutoff in 1.5 2 3; do
			Rscript $writescript $dist $cutoff $method $model&
				if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
					# define maxjobs and n using maxjobsn skeleton
				    wait # wait until all have finished (not optimal, but most times good enough)
				    echo $n wait
				fi
		done
	done
done
wait
