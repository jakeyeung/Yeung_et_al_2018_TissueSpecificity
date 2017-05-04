#!/bin/sh
# Jake Yeung
# run.run_mara_on_liver_specific_DHS_SqlDiffModsDiffCutoffs.sh
# Run run_mara_on_liver_specific_DHS_SqlDiffModsDiffCutoffs.R with different parameters 
# 2017-04-29

# rscript="/home/yeung/projects/tissue-specificity/scripts/three_way_contingency_tables/run_mara_on_liver_specific_DHS_SqlDiffModsDiffCutoffs.R"
# rscript="/home/yeung/projects/tissue-specificity/scripts/three_way_contingency_tables/run_mara_on_liver_specific_DHS_SqlDiffModsDiffCutoffs.byargs.R"
# sumscript="/home/yeung/projects/tissue-specificity/scripts/three_way_contingency_tables/downstream.run_mara_on_liver_specific_DHS.byargs.R"
runscript="/home/yeung/projects/tissue-specificity/scripts/three_way_contingency_tables/run.mara_and_downstream.byargs.sh"
[[ ! -e $runscript ]] && echo "$runscript not found, exiting" && exit 1

distfilt=40000
jcutoff="2"
jcutoff="1.5"
cl="0.5"
jweight=0
inclproms="FALSE"
dopairs="FALSE"

# jmod="Liver_SV129"
# jcutofflow="0"

for jmod in "Liver_SV129" "Liver_SV129,Liver_BmalKO" "Kidney_SV129" "Kidney_SV129,Kidney_BmalKO"; do
	bash $runscript $jweight $distfilt $jmod $jcutoff $cl $inclproms $dopairs&
	# echo $cl
	# Rscript $rscript $jweight $distfilt $jmod $jcutoff $cl $inclproms
	# ret=$?; [[ $ret -ne 0  ]] && echo "ERROR: script failed" && exit 1
	# Rscript $sumscript $jweight $distfilt $jmod $jcutoff $cl $inclproms
	# ret=$?; [[ $ret -ne 0  ]] && echo "ERROR: script failed" && exit 1
done
wait
# ret=$?; [[ $ret -ne 0  ]] && echo "ERROR: script failed" && exit 1

