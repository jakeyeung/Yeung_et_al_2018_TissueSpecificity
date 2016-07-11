#!/bin/sh
# Jake Yeung
# run.write_gene_exprs_table_for_mara.modules.sh
# Run write table in batch
# 2016-07-10


writescript="/home/yeung/projects/tissue-specificity/scripts/liver_kidney_WTKO/write_gene_exprs_table_for_mara.modules.R"
[[ ! -e $writescript ]] && echo "$writescript not found, exiting" && exit 1
# need: jmodel, outdir, jmethod

for jmethod in "BIC" "g=1001"; do
	for jmodel in "Liver_SV129" "Kidney_SV129" "Liver_SV129,Kidney_SV129"; do
		outdir="/home/yeung/projects/tissue-specificity/data/gene_exprs/liver_v_kidney/atger_with_kidney.$jmodel.$jmethod"
		Rscript $writescript $jmodel $outdir $jmethod
	done
done
