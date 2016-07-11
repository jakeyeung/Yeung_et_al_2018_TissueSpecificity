#!/bin/sh
# Jake Yeung
# copy_genelist_filter_motifs_copy_motifs_from_vitalit.sh
# Combine several scripts together to copy genelist to vitalit, filter motifs by genelist. copy motifs back into Rstdudio 
# 2016-06-27

remotehost=jyeung@frt.el.vital-it.ch

glistscript="/home/yeung/projects/tissue-specificity/shellscripts/copy_genelist_filter_motifs_copy_motifs/util.copy_genelist_to_vitalit.sh"
filtscriptVITAL="/Home/jyeung/projects/tissue_specificity_hogenesch_shellscripts/motevo_dhs_scripts_clean/3-run_motevo/6d-filter_by_gene_list.sh"
motifscript="/home/yeung/projects/tissue-specificity/shellscripts/copy_genelist_filter_motifs_copy_motifs/util.copy_motifs_under_peaks_from_vitalit.sh"

[[ ! -e $glistscript ]] && echo "$glistscript not found, exiting" && exit 1
[[ ! -e $motifscript ]] && echo "$motifscript not found, exiting" && exit 1

# copy gene list to vitalit
genelistdir="/home/yeung/projects/tissue-specificity/data/gene_lists/liver_kidney_wtko_modules_allgenes"
genelistoutfVITAL="/scratch/el/monthly/jyeung/motevo_dhs_outputs/gene_lists/liver_kidney_wtko/liver_kidney_wtmodules.all_genes.bic_g1001.txt"
bash $glistscript $genelistdir $genelistoutfVITAL
ret=$?; [[ $ret -ne 0  ]] && echo "ERROR: script failed" && exit 1

# filter script on vitalit
filtdirVITAL="/scratch/el/monthly/jyeung/motevo_dhs_outputs/motevo_outputs/liver_kidney_wtko_closestbed_multiple_genes_ncond_modules.all_genes.bic_g1001"
ssh $remotehost "bash $filtscriptVITAL $genelistoutfVITAL $filtdirVITAL"
ret=$?; [[ $ret -ne 0  ]] && echo "ERROR: script failed" && exit 1

# copy bed files into rstudio
outfbed="/home/yeung/data/tissue_specificity/motevo_dhs/liver_kidney_wtko_modules/all_sites.closest.filter.genelst.dist.50000.all_genes.bic_g1001.bed"
bash $motifscript $filtdirVITAL $outfbed
ret=$?; [[ $ret -ne 0  ]] && echo "ERROR: script failed" && exit 1

