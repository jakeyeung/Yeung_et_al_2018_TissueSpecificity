# 2016-08-28
# Jake Yeung
# Summarize differences in phases between mRNA and regulator. Include an alpha value?

rm(list=ls())

library(ggplot2)
library(hash)

jmeth <- "g=1001"

source("scripts/functions/ProteomicsFunctions.R")
source("scripts/functions/GetTFs.R")
source("scripts/functions/PlotFunctions.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/PlotActivitiesFunctions.R")
source("scripts/functions/LoadActivitiesLong.R")
source("scripts/functions/SvdFunctions.R")

MatchMotifToGene <- function(motif, tfs, genes.lst.filt, fits.bytiss.sub){
  # Match motif to genes, given a gene list to filter (to handle multiple genes)
  genes.all <- GetGenesFromMotifs(motif, tfs)
  genes.filt <- genes.all[which(genes.all %in% genes.lst.filt)]
  if (length(genes.filt) == 0){
    jgene <- NA
  } else if (length(genes.filt) == 1){
    jgene <- genes.filt[[1]]
  } else {
    # more than one gene, pick gene with better pvalue
    fits.sub <- subset(fits.bytiss.sub, gene %in% genes.filt)
    jgene <- as.character(fits.sub$gene[which(fits.sub$pval == min(fits.sub$pval))])
    warning(paste("Multiple genes, picking best fit", jgene))
  }
  return(jgene)
}

MergeActivityProteinRNA <- function(jmod, fits.long.filt, fits.bytiss, prot, jtiss = "Liver_SV129"){
  genes.liv <- as.character(subset(fits.long.filt, model == jmod)$gene)
  
  tfs <- GetTFs(split.commas = TRUE, get.mat.only = FALSE)
  tfs.mat <- GetTFs(split.commas = TRUE, get.mat.only = TRUE)
  genes.tfs <- genes.liv[genes.liv %in% tfs]
  # genes.tfs <- genes.liv
  
  # Can we get TFs that do not have motifs? ---------------------------------
  # see if Jingkui have already annotated proteins as TF or not. 
  
  fits.sub <- subset(fits.bytiss, gene %in% genes.tfs & tissue == jtiss, select = c(gene, amp, phase, pval, exper))
  prot.sub <- subset(prot, gene %in% genes.tfs, select = c(gene, amp, phase, pval, exper))
  # if gene assigned to two peptides, take one with best rhythmicity
  prot.sub <- prot.sub %>%
    group_by(gene) %>%
    filter(pval == min(pval))
  merged.long <- rbind(as.data.frame(fits.sub), as.data.frame(prot.sub))
  merged <- dcast(merged.long, formula = gene ~ exper, value.var = "phase")
  # add worst pvalue
  merged.worstpval <- merged.long %>%
    group_by(gene) %>%
    filter(pval == max(pval))
  worstpvals <- hash(as.character(merged.worstpval$gene), merged.worstpval$pval)
  
  # add to merged
  merged$pval.worst <- sapply(as.character(merged$gene), function(g) worstpvals[[g]])
  
  outbase <- "/home/yeung/projects/tissue-specificity/results/MARA.liver_kidney"
  outmain <- file.path(outbase, paste0("promoters.", jmod, ".g=1001"))
  indir <- file.path(outmain, "atger_with_kidney.bugfixed")
  act.long <- LoadActivitiesLongKidneyLiver(indir, collapse.geno.tissue = TRUE, shorten.motif.name = TRUE)
  
  omega <- 2 * pi / 24
  act.complex <- subset(act.long, tissue == jtiss) %>%
    group_by(gene, tissue) %>%
    do(ProjectToFrequency2(., omega, add.tissue=TRUE))
  act.complex$phase <- ConvertArgToPhase(Arg(act.complex$exprs.transformed), omega)
  
  act.complex$motif <- act.complex$gene
  act.complex$gene <- sapply(as.character(act.complex$motif), MatchMotifToGene, tfs.mat, as.character(merged.long$gene), subset(fits.bytiss, tissue == jtiss))
  act.sub <- subset(act.complex, !is.na(gene))
  phase.hash <- hash(act.sub$gene, act.sub$phase)
  # annotate merged
  print(merged)
  merged$activity <- sapply(as.character(merged$gene), function(g){
    act.phase <- phase.hash[[g]]
    if (is.null(act.phase)){
      act.phase <- NA
    }
    return(act.phase)
  })
  return(merged)
}

source("scripts/functions/GetClockGenes.R")


# Load --------------------------------------------------------------------

load("Robjs/liver_kidney_atger_nestle/dat.long.liverkidneyWTKO.bugfixed.Robj", v=T)
load("Robjs/liver_kidney_atger_nestle/fits.long.multimethod.filtbest.staggeredtimepts.bugfixed.annotated.Robj", v=T)
load("Robjs/liver_kidney_atger_nestle/fits.bytiss.bugfixed.Robj", v=T)
fits.long.filt <- subset(fits.long.filt, method == jmeth)
prot.long <- LoadProteomicsData()

prot <- LoadProteomicsData(as.long = FALSE)
prot <- subset(prot, select = c(Gene.names, mean, amp, relamp, phase, pval, qv))
colnames(prot)[colnames(prot) == "Gene.names"] <- "gene"
prot$exper <- "prot"
fits.bytiss$exper <- "rnaseq"

merged.test <- MergeActivityProteinRNA(jmod = "Liver_SV129,Liver_BmalKO", fits.long.filt, fits.bytiss, prot, jtiss = "Liver_SV129")

# Get rhythmic TFs and their motifs ---------------------------------------

jmod <- "Liver_SV129,Liver_BmalKO"
genes.liv <- as.character(subset(fits.long.filt, model == jmod)$gene)

tfs <- GetTFs(split.commas = TRUE, get.mat.only = FALSE)
tfs.mat <- GetTFs(split.commas = TRUE, get.mat.only = TRUE)
genes.tfs <- genes.liv[genes.liv %in% tfs]
# genes.tfs <- genes.liv

# Can we get TFs that do not have motifs? ---------------------------------
# see if Jingkui have already annotated proteins as TF or not. 

fits.sub <- subset(fits.bytiss, gene %in% genes.tfs & tissue == "Liver_SV129", select = c(gene, amp, phase, pval, exper))
prot.sub <- subset(prot, gene %in% genes.tfs, select = c(gene, amp, phase, pval, exper))
# if gene assigned to two peptides, take one with best rhythmicity
prot.sub <- prot.sub %>%
  group_by(gene) %>%
  filter(pval == min(pval))
merged.long <- rbind(as.data.frame(fits.sub), as.data.frame(prot.sub))
merged <- dcast(merged.long, formula = gene ~ exper, value.var = "phase")
# add worst pvalue
merged.worstpval <- merged.long %>%
  group_by(gene) %>%
  filter(pval == max(pval))
worstpvals <- hash(as.character(merged.worstpval$gene), merged.worstpval$pval)

# add to merged
merged$pval.worst <- sapply(as.character(merged$gene), function(g) worstpvals[[g]])
merged.sub <- subset(merged, pval.worst < 0.05)

ggplot(merged.sub, aes(x = rnaseq, y = prot, label = gene)) +
  geom_text() + theme_bw(18) + geom_abline(slope = 1) + theme(legend.position = "bottom")

ggplot(merged, aes(x = rnaseq, y = prot, label = gene, alpha = -log10(pval.worst))) +
  geom_text() + theme_bw(18) + geom_abline(slope = 1) + theme(legend.position = "bottom")
  

# Add Activity ------------------------------------------------------------

jmod <- "Liver_SV129,Liver_BmalKO"
outbase <- "/home/yeung/projects/tissue-specificity/results/MARA.liver_kidney"
outmain <- file.path(outbase, paste0("promoters.", jmod, ".g=1001"))
indir <- file.path(outmain, "atger_with_kidney.bugfixed")

act.long <- LoadActivitiesLongKidneyLiver(indir, collapse.geno.tissue = TRUE, shorten.motif.name = TRUE)

omega <- 2 * pi / 24
act.complex <- act.long %>%
  group_by(gene, tissue) %>%
  do(ProjectToFrequency2(., omega, add.tissue=TRUE))
act.complex$phase <- ConvertArgToPhase(Arg(act.complex$exprs.transformed), omega)


act.complex$motif <- act.complex$gene
genes <- sapply(as.character(act.complex$motif), MatchMotifToGene, tfs.mat, as.character(merged.long$gene))
act.complex$gene <- genes

s.act <- SvdOnComplex(act.complex, value.var = "exprs.transformed")


jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)

# plot motifs
max.labs <- 20
jtitle <- ""
comp <- 1
eigens.act <- GetEigens(s.act, period = 24, comp = comp, adj.mag = TRUE, constant.amp = 4, label.n = max.labs, jtitle = jtitle, peak.to.trough = TRUE, label.gene = c("bHLH_family.p2", "RORA.p2", "SRF.p3", "HSF1.2.p2", "TFAP2B.p2"))
print(eigens.act$u.plot + ggtitle(jmod))

PlotActivitiesWithSE(subset(act.long, gene == "CEBPA.B_DDIT3.p2"))
PlotActivitiesWithSE(subset(act.long, gene == "POU2F1..3.p2"))
PlotActivitiesWithSE(subset(act.long, gene == "RFX1..5_RFXANK_RFXAP.p2"))
PlotActivitiesWithSE(subset(act.long, gene == "HNF4A_NR2F1.2.p2"))

PlotmRNAActivityProtein(dat.long, act.long, prot.long, "Mafb", "MAFB.p2", "Mafb")
PlotmRNAActivityProtein(dat.long, act.long, prot.long, "Cebpb", "CEBPA.B_DDIT3.p2", "Cebpb")
PlotmRNAActivityProtein(dat.long, act.long, prot.long, "Egr1", "EGR1..3.p2", "Egr1")
PlotmRNAActivityProtein(dat.long, act.long, prot.long, "Pou2f1", "POU2F1..3.p2", "Pou2f1")
PlotmRNAActivityProtein(dat.long, act.long, prot.long, "Rfxap", "RFX1..5_RFXANK_RFXAP.p2", "Rfxap")
PlotmRNAActivityProtein(dat.long, act.long, prot.long, "Hnf4a", "HNF4A_NR2F1.2.p2", "Hnf4a")
PlotmRNAActivityProtein(dat.long, act.long, prot.long, "Jun", "JUN.p2", "Jun")
PlotmRNAActivityProtein(dat.long, act.long, prot.long, "Creb3", "ATF5_CREB3.p2", "Creb3")
