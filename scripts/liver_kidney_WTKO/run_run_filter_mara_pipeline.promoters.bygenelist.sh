#!/bin/sh
# Jake Yeung
# run_run_filter_mara_pipeline.promoters.bygenelist.sh
# Using genelist
# 2016-04-12

initscript="/home/yeung/projects/tissue-specificity/scripts/liver_kidney_WTKO/write_gene_exprs_table_for_mara.modules.R"
combinescript="/home/yeung/projects/ridge-regression/run_scripts/combine_activities_and_plot.one_dir.sh"
runscript="/home/yeung/projects/ridge-regression/run_scripts/run_mara_batch_promoters.sh"

# geneexprsmain="/home/yeung/projects/tissue-specificity/data/gene_exprs/liver_v_kidney"
geneexprsdir="/home/yeung/projects/tissue-specificity/data/gene_exprs/liver_v_kidney/atger_with_kidney"
outmain="/home/yeung/projects/tissue-specificity/results/MARA.liver_kidney"
[[ ! -d $outmain ]] && echo "$outmain not found, exiting" && exit 1

# Nmain="/home/yeung/projects/tissue-specificity/data/sitecounts/motevo/liver_kidney_sitecounts"
# Nmain="/home/yeung/projects/tissue-specificity/data/sitecounts/motevo/liver_kidney_sitecounts_tissuespecpeaks"
Nmain="/home/yeung/projects/tissue-specificity/data/sitecounts/motevo/liver_kidney_sitecounts_tissuespecpeaks_cutofflow"
# dist=40000
cutofflow="0.5"
cutoff=3
cross="TRUE"  # or FALSE
tspeaks="TRUE"
# method="BIC"
# for model in "Liver_SV129" "Kidney_SV129" "Liver_SV129,Kidney_SV129"; do
for cutoff in "1.5" "2" "3"; do
	for dist in 10000 20000 40000; do
		for method in "BIC" "g=1001"; do
			for model in "Liver_SV129"; do
				echo $model
				# outdir="/home/yeung/projects/tissue-specificity/results/MARA/atger_kidney_enhancers_40kb_cutoff3.$model"
				outdir=$outmain/atger_kidney.dist_"$dist".cutoff_"$cutoff".model_"$model".method_"$method".cross_"$cross".tspeaks_"$tspeaks".cutofflow_"$cutofflow"

				# Gene expression directory should be all genes if we use an N matrix that is filtered for genes
				# geneexprsdir="/home/yeung/projects/tissue-specificity/data/gene_exprs/liver_v_kidney/atger_with_kidney.$model"
				# geneexprsdir="$geneexprsmain/atger_kidney.$model.$method"
				subdir=$(basename $geneexprsdir)

				# Nmat="/home/yeung/projects/tissue-specificity/data/sitecounts/motevo/sitecount_matrix_geneids"
				# Nmat="/home/yeung/projects/tissue-specificity/data/sitecounts/motevo/liver_kidney_sitecounts/liver_rhythmic_genes_methodbic_cutoff3_distance40000.mat"
				# Nmat="/home/yeung/projects/tissue-specificity/data/sitecounts/motevo/liver_kidney_sitecounts/liver_rhythmic_genes.method_BIC.dist_40000.cutoff_3.mat"
				Nmat=$Nmain/sitecounts_enhancers.model_"$model".method_"$method".dist_"$dist".cutoff_"$cutoff".cross_"$cross".tspeaks_"$tspeaks".cutofflow_"$cutofflow".mat

				[[ -d $outdir ]] && echo "$outdir found, exiting" && continue  # overwriting mode does not work! need this to prevent bugs: TODO just remove safely the directory?
				[[ ! -d $outdir ]] && mkdir $outdir
				[[ ! -e $runscript ]] && echo "$runscript not found, exiting" && exit 1
				[[ ! -e $Nmat ]] && echo "$Nmat not found, exiting" && exit 1

				# generate gene expression directory

				# do not need to generate this, we have it already
				# echo "Rscript $initscript $model $geneexprsdir $method"
				# Rscript $initscript $model $geneexprsdir $method
				# ret=$?; [[ $ret -ne 0  ]] && echo "ERROR: script failed" && exit 1
				[[ ! -d $geneexprsdir ]] && echo "$geneexprsdir not found, exiting" && exit 1

				# run MARA here
				bash $runscript $geneexprsdir $outdir $Nmat
				ret=$?; [[ $ret -ne 0  ]] && echo "ERROR: script failed" && exit 1
				bash $combinescript $outdir/$subdir
				ret=$?; [[ $ret -ne 0  ]] && echo "ERROR: script failed" && exit 1
			done
		done
	done
done

