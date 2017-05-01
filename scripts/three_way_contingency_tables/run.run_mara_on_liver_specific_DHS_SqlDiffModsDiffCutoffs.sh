#!/bin/sh
# Jake Yeung
# run.run_mara_on_liver_specific_DHS_SqlDiffModsDiffCutoffs.sh
# Run run_mara_on_liver_specific_DHS_SqlDiffModsDiffCutoffs.R with different parameters 
# 2017-04-29

# rscript="/home/yeung/projects/tissue-specificity/scripts/three_way_contingency_tables/run_mara_on_liver_specific_DHS_SqlDiffModsDiffCutoffs.R"
rscript="/home/yeung/projects/tissue-specificity/scripts/three_way_contingency_tables/run_mara_on_liver_specific_DHS_SqlDiffModsDiffCutoffs.byargs.R"
distfilt=40000
jcutoff="2"
# jmod="Liver_SV129"
# jcutofflow="0"

jweight=0
# for jmod in "Kidney_SV129" "Kidney_SV129;Kidney_BmalKO"; do
for jmod in "Liver_SV129" "Liver_SV129,Liver_BmalKO"; do
	for cl in 0 0.5 1.5 2; do
		echo $cl
		Rscript $rscript $jweight $distfilt $jmod $jcutoff $cl
		ret=$?; [[ $ret -ne 0  ]] && echo "ERROR: script failed" && exit 1
	done
done
