#!/bin/sh
# Jake Yeung
# run.run_plda_with_sql_db.sh
# Run plda from sql db in order to query large database quickly rather than generating the Ns 
# 2016-09-14

rscript="/home/yeung/projects/tissue-specificity/scripts/sql_db_scripts/run_plda_with_sql_db.R"

for m in "Liver_SV129" "Liver_SV129,Liver_BmalKO" "Liver_BmalKO" "Kidney_SV129" "Kidney_SV129,Kidney_BmalKO"; do
# for m in "Liver_SV129,Liver_BmalKO"; do
# for m in "Liver_SV129"; do
	Rscript $rscript $m
done
