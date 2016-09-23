#!/bin/sh
# Jake Yeung
# run.calculate_liver_specific_peaks_by_genelists.sh
# Run Rscript
# 2016-08-09

jscript="/home/yeung/projects/tissue-specificity/scripts/liver_kidney_WTKO/calculate_liver_specific_peaks_by_genelists.R"

jtiss="Liver"
for mod in "Liver_SV129" "Liver_SV129,Liver_BmalKO" "Kidney_SV129" "Kidney_SV129,Kidney_BmalKO"; do
	echo $mod
	Rscript $jscript $mod $jtiss 14&
	# echo "Rscript $jscript $mod 1000"
done
wait

# jtiss="Kidney"
# for mod in "Kidney_SV129" "Kidney_SV129,Kidney_BmalKO" "Liver_SV129" "Liver_SV129,Liver_BmalKO"; do
# 	echo $mod
# 	Rscript $jscript $mod $jtiss 15&
# 	# echo "Rscript $jscript $mod 1000"
# done
# wait
