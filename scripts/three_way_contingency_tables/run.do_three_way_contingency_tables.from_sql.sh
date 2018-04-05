#!/bin/sh
# Jake Yeung
# run.do_three_way_contingency_tables.from_sql.sh
# Run it  
# 2017-05-08

rscript="/home/yeung/projects/tissue-specificity/scripts/three_way_contingency_tables/do_three_way_contingency_tables.from_sql.R"

for nblevel in 2 3; do
	for cutoff in 0 0.8; do
		Rscript $rscript $nblevel $cutoff
	done
done
