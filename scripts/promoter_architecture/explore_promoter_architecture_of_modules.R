# 2016-08-15
# Jake Yeung
# Investigate whether CpG and TATA can separate between modules

rm(list=ls())

library(dplyr)
library(ggplot2)
library(hash)

source('scripts/functions/PlotGeneAcrossTissues.R')

# Load S1 from westermark paper -------------------------------------------

inf <- "/home/shared/promoters/promoters.westermark_s1_table/westermark_s1_table.fixed.txt"

annots <- read.table(inf, sep = "\t", header = TRUE)


# Load modules ------------------------------------------------------------

load("Robjs/dat.long.fixed_rik_genes.Robj", v=T); dat.hog <- dat.long; rm(dat.long)
load("Robjs/liver_kidney_atger_nestle/dat.long.liverkidneyWTKO.bugfixed.Robj", v=T)
load("Robjs/liver_kidney_atger_nestle/fits.long.multimethod.filtbest.staggeredtimepts.bugfixed.annotated.Robj", v=T)
load("Robjs/nconds_g1000_11_tissues/fits_long.11_tiss_3_max.g1000.bestmodel.filteramp.0.15.Robj", v=T)
fits.long.filt <- subset(fits.long.filt, method == "g=1001")


# Compare liver modules ---------------------------------------------------

livWT.genes <- as.character(subset(fits.long.filt, model == "Liver_SV129")$gene)
livWTKO.genes <- as.character(subset(fits.long.filt, model == "Liver_SV129,Liver_BmalKO")$gene)
flat.genes <- as.character(subset(fits.long.filt, model == "")$gene)
tw.genes <- as.character(subset(fits.long, n.rhyth >= 8)$gene)

# check overlap of tw.genes with liverWT and liverWTKO genes.
tw.livWT.overlaps <- intersect(tw.genes, livWT.genes)
tw.livWTKO.overlaps <- intersect(tw.genes, livWTKO.genes)

tw.genes.nonoverlap <- tw.genes[! tw.genes %in% tw.livWT.overlaps & ! tw.genes %in% tw.livWTKO.overlaps]

# annotate annots table with livWT vs livWTKO genes
mods.hash <- hash(as.character(fits.long.filt$gene), as.character(fits.long.filt$model))
is.tw.hash <- hash(tw.genes.nonoverlap, TRUE)

jmods <- c("Flat", "Liver_SV129", "Liver_SV129,Liver_BmalKO", "Kidney_SV129")

annots$model <- sapply(as.character(annots$Gene.symbol), function(g){
  if (g == ""){
    g <- "BlankGene"
  }
  jmodel <- mods.hash[[g]]
  if (is.null(jmodel)){
    jmodel <- NA
  } else if (jmodel == ""){
    jmodel <- "Flat"
  }
  if (! jmodel %in% jmods){
    # try to assign to TW or not
    is.tw <- is.tw.hash[[g]]
    if (!is.null(is.tw)){
      jmodel <- "TissueWide"
    }
  }
  return(jmodel)
}, USE.NAMES = TRUE)



# Compare CpG ratios ------------------------------------------------------

jmods.tw <- c(jmods, "TissueWide")
ggplot(subset(annots, model %in% jmods.tw), aes(x = model, y = CpG.ratio)) + geom_boxplot() 
ggplot(subset(annots, model %in% jmods.tw), aes(x = model, y = PI)) + geom_boxplot() 
ggplot(subset(annots, model %in% jmods.tw), aes(x = model, y = Half.life..hr.)) + geom_boxplot() 
ggplot(subset(annots, model %in% jmods.tw), aes(x = model, y = H3K4me3..TSS)) + geom_boxplot() 
ggplot(subset(annots, model %in% jmods.tw), aes(x = model, y = Nucleosome..1.H2A.Z)) + geom_boxplot() 

# summary of TATA boxes
annots.sum <- subset(annots, model %in% jmods.tw) %>%
  group_by(model) %>%
  summarise(n.TATAbox = length(which(TATA.box == TRUE)),
            n.BREu = length(which(BREu == TRUE)),
            n.BREd = length(which(BREd == TRUE)),
            n.CTFbinding = length(which(CTF.binding == TRUE)),
            n.genes = length(Gene.symbol)) %>%
  mutate(frac.TATAbox = n.TATAbox/n.genes,
         frac.BREu = n.BREu/n.genes,
         frac.BREd = n.BREd/n.genes,
         frac.CTFbinding = n.CTFbinding/n.genes)

ggplot(annots.sum, aes(x = model, y = frac.TATAbox)) + geom_bar(stat = "identity")
ggplot(annots.sum, aes(x = model, y = frac.BREu)) + geom_bar(stat = "identity")
ggplot(annots.sum, aes(x = model, y = frac.BREd)) + geom_bar(stat = "identity")
ggplot(annots.sum, aes(x = model, y = frac.CTFbinding)) + geom_bar(stat = "identity")
