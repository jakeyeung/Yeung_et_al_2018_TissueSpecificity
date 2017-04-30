#!/bin/sh
# Jake Yeung
# run.downstream.run_mara_on_liver_specific_DHS.byargs.sh
#  
# 2017-04-29

sumscript="/home/yeung/projects/tissue-specificity/scripts/three_way_contingency_tables/downstream.run_mara_on_liver_specific_DHS.byargs.R"

[[ ! -e $sumscript ]] && echo "$sumscript not found, exiting" && exit 1

jweight=0
distfilt=40000
# jmod="Liver_SV129"
jcutoff="2"

# for jmod in "Kidney_SV129"; do
# for jcutofflow in 0 0.5 1.5 2; do
for jmod in "Liver_SV129" "Liver_SV129,Liver_BmalKO"; do
	for jcutoff in "2" "3"; do
		for jcutofflow in 0 0.5 1.5 2 2.5 3; do
			echo $jcutofflow
			Rscript $sumscript $jweight $distfilt $jmod $jcutoff $jcutofflow
			ret=$?; [[ $ret -ne 0  ]] && echo "ERROR: script failed" && exit 1
		done
	done
done

