#!/bin/sh
# Jake Yeung
# write_sitecounts_run_mara_pipeline.sh
# Gluing write N sitecounts and run MARA for many parameters 
# 2016-07-12

Nmain="/home/yeung/projects/tissue-specificity/data/sitecounts/motevo/liver_kidney_sitecounts_tissuespecpeaks_cutofflow0_singletons"
maraout="/home/yeung/projects/tissue-specificity/results/MARA.liver_kidney/singletons"
plotdir="/home/yeung/projects/tissue-specificity/plots/MARA/singletons_tspeak"
plotout=$plotdir/"singletons.tspeak.pdf"


# write sitecounts
[[ ! -d $Nmain ]] && mkdir $Nmain
writescript="/home/yeung/projects/tissue-specificity/scripts/liver_kidney_WTKO/run.write_N_sitecounts_table_for_mara.sh"
echo "Rscript $writescript $Nmain"
bash $writescript $Nmain
ret=$?; [[ $ret -ne 0  ]] && echo "ERROR: script failed" && exit 1
wait

# run MARA
[[ ! -d $maraout ]] && mkdir $maraout
marascript="/home/yeung/projects/tissue-specificity/scripts/liver_kidney_WTKO/run_run_filter_mara_pipeline.promoters.bygenelist.sh"
bash $marascript $Nmain $maraout
# ret=$?; [[ $ret -ne 0  ]] && echo "ERROR: script failed" && exit 1
wait

# plot output
[[ ! -d $plotdir ]] && mkdir $plotdir
plotscript="/home/yeung/projects/tissue-specificity/scripts/liver_kidney_WTKO/run.analyze_activities.batch.R"
Rscript $plotscript $maraout $plotout
ret=$?; [[ $ret -ne 0  ]] && echo "ERROR: script failed" && exit 1
wait

