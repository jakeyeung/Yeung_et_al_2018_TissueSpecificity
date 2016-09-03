#!/bin/sh
# Jake Yeung
# run_run_filter_mara_pipeline.promoters.bygenelist.sh
# Using genelist
# 2016-04-12

# Nmain="/home/yeung/projects/tissue-specificity/data/sitecounts/motevo/liver_kidney_sitecounts_tissuespecpeaks_cutofflow_limited_crossprods"
# outmain="/home/yeung/projects/tissue-specificity/results/MARA.liver_kidney.limitprod"
Nmat=$1
outdir=$2
geneexprsdir=$3
lambda=$4
echo "Nmat: $Nmat"
echo "Outdir: $outdir"
echo "Gene exprs: $geneexprsdir"
echo "Lambda: $lambda"

# if not supplied then assume liver and kidney
if [ -z $geneexprsdir ]; then
	geneexprsdir="/home/yeung/projects/tissue-specificity/data/gene_exprs/liver_v_kidney/atger_with_kidney.bugfixed"
fi

initscript="/home/yeung/projects/tissue-specificity/scripts/liver_kidney_WTKO/write_gene_exprs_table_for_mara.modules.R"
combinescript="/home/yeung/projects/ridge-regression/run_scripts/combine_activities_and_plot.one_dir.sh"
runscript="/home/yeung/projects/ridge-regression/run_scripts/run_mara_batch_promoters.sh"
[[ ! -e $runscript ]] && echo "Runscript $runscript not found, exiting" && exit 1

# geneexprsmain="/home/yeung/projects/tissue-specificity/data/gene_exprs/liver_v_kidney"
# geneexprsdir="/home/yeung/projects/tissue-specificity/data/gene_exprs/liver_v_kidney/atger_with_kidney.bugfixed"

bname=$(basename $Nmat)
subdir=$(basename $geneexprsdir)
[[ -d $outdir ]] && echo "Outdir $outdir found, exiting" && exit 1  # overwriting mode does not work! need this to prevent bugs: TODO just remove safely the directory?
[[ ! -d $outdir ]] && mkdir $outdir
[[ ! -e $Nmat ]] && echo "Nmat $Nmat not found, exiting" && exit 1
[[ ! -d $geneexprsdir ]] && echo "Geneexprsdir $geneexprsdir not found, exiting" && exit 1
# run MARA here
bash $runscript $geneexprsdir $outdir $Nmat $lambda
ret=$?; [[ $ret -ne 0  ]] && echo "ERROR: script failed" && exit 1
bash $combinescript $outdir/$subdir
ret=$?; [[ $ret -ne 0  ]] && echo "ERROR: script failed" && exit 1
