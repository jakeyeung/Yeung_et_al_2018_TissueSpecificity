#!/bin/sh
# Jake Yeung
# run.mara_and_downstream.byargs.sh
# Run MARA and Downstream in a single script  
# 2017-05-02


rscript="/home/yeung/projects/tissue-specificity/scripts/three_way_contingency_tables/run_mara_on_liver_specific_DHS_SqlDiffModsDiffCutoffs.byargs.R"
sumscript="/home/yeung/projects/tissue-specificity/scripts/three_way_contingency_tables/downstream.run_mara_on_liver_specific_DHS.byargs.R"
[[ ! -e $sumscript ]] && echo "$sumscript not found, exiting" && exit 1
[[ ! -e $rscript ]] && echo "$rscript not found, exiting" && exit 1

# distfilt=40000
# jcutoff="1.75"
# jcutoff="2"
# # jmod="Liver_SV129"
# # jcutofflow="0"
# 
# jweight=0
# inclproms="TRUE"
# 
# jcutofflow="0.5"

jweight=$1
distfilt=$2
jmod=$3
jcutoff=$4
jcutofflow=$5
inclproms=$6

Rscript $rscript $jweight $distfilt $jmod $jcutoff $jcutofflow $inclproms
ret=$?; [[ $ret -ne 0  ]] && echo "ERROR: script failed" && exit 1
Rscript $sumscript $jweight $distfilt $jmod $jcutoff $jcutofflow $inclproms
ret=$?; [[ $ret -ne 0  ]] && echo "ERROR: script failed" && exit 1
