# 2016-08=08
# how many liver-specific peaks there are in different sets? Redo with liver kidney WT KO 

rm(list=ls())

library(ggplot2)

setwd("/home/yeung/projects/tissue-specificity")

library(dplyr)
library(ggplot2)
library(parallel)
library(hash)

source("scripts/functions/PlotGeneAcrossTissues.R")

load("Robjs/liver_kidney_atger_nestle/fits.long.multimethod.filtbest.staggeredtimepts.bugfixed.Robj", v=T)
load("Robjs/liver_kidney_atger_nestle/dat.long.liverkidneyWTKO.bugfixed.Robj", v=T)
load("Robjs/liver_kidney_atger_nestle/fits.bytiss.bugfixed.Robj", v=T)
fits.long.filt <- subset(fits.long.filt, method == "g=1001")


# Load output from calculate_liver_specific_peaks -------------------------

# load("/home/yeung/projects/tissue-specificity/Robjs/liver_kidney_atger_nestle/tissue_specific_peaks/n.tiss.spec.df.out.lst.rand.100.tissue.Liver.module.Liver_SV129.Robj", v=T)
# wtmodulef <- "/home/yeung/projects/tissue-specificity/Robjs/liver_kidney_atger_nestle/tissue_specific_peaks/n.tiss.spec.df.out.lst.rand.100.tissue.Liver.module.Liver_SV129random.flat.TRUE.Robj"
wtmodulef <- "/home/yeung/projects/tissue-specificity/Robjs/liver_kidney_atger_nestle/tissue_specific_peaks/n.tiss.spec.df.out.lst.rand.1000.tissue.Liver.module.Liver_SV129random.flat.TRUE.Robj"
load(wtmodulef, v=T)
df.out.lst.merged.liverWT <- df.out.lst.merged
df.out.lst.merged.liverWT$gene.type[1:3] <- c("Liver_SV129", "Flat_SV129", "Flat.filt_SV129")

# wtko.modulef <- "/home/yeung/projects/tissue-specificity/Robjs/liver_kidney_atger_nestle/tissue_specific_peaks/n.tiss.spec.df.out.lst.rand.100.tissue.Liver.module.Liver_SV129,Liver_BmalKO.Robj"
# wtko.modulef <- "/home/yeung/projects/tissue-specificity/Robjs/liver_kidney_atger_nestle/tissue_specific_peaks/n.tiss.spec.df.out.lst.rand.100.tissue.Liver.module.Liver_SV129,Liver_BmalKOrandom.flat.TRUE.Robj"
wtko.modulef <- "/home/yeung/projects/tissue-specificity/Robjs/liver_kidney_atger_nestle/tissue_specific_peaks/n.tiss.spec.df.out.lst.rand.1000.tissue.Liver.module.Liver_SV129,Liver_BmalKOrandom.flat.TRUE.Robj"
load(wtko.modulef, v=T)
df.out.lst.merged.liverWTKO <- df.out.lst.merged
df.out.lst.merged.liverWTKO$gene.type[1:3] <- c("Liver_SV129.Liver_BmalKO", "Flat_SV129BmalKO", "Flat.filt_SV129BmalKO")

df.out.lst.merged <- rbind(df.out.lst.merged.liverWT, df.out.lst.merged.liverWTKO)

df.out.lst.merged$xlabs <- make.names(df.out.lst.merged$gene.type, unique = TRUE)

df.out.lst.merged <- df.out.lst.merged[order(df.out.lst.merged$gene.type, decreasing = FALSE), ]

# Show barplot? -----------------------------------------------------------

ggplot(df.out.lst.merged, aes(x = xlabs, y = mean.tiss.spec.per.gene)) + geom_bar(stat = "identity")
ggplot(df.out.lst.merged, aes(x = xlabs, y = frac.n.spec.by.gene)) + geom_bar(stat = "identity")

# Show significance  ------------------------------------------------------

df.out.lst.bg <- df.out.lst.merged[grepl("^Random", df.out.lst.merged$xlabs), ]
df.out.lst.bg.meanvar <- df.out.lst.bg %>%
  group_by(gene.type) %>%
  summarise(mean.frac = mean(frac.n.spec.by.gene),
            sd.frac = sd(frac.n.spec.by.gene))
df.out.lst.fg.meanvar <- df.out.lst.merged[!grepl("^Random", df.out.lst.merged$xlabs), ] %>%
  group_by(gene.type) %>%
  summarise(mean.frac = mean(frac.n.spec.by.gene),
            sd.frac = sd(frac.n.spec.by.gene))
df.out.lst.meanvar <- rbind(df.out.lst.fg.meanvar[!grepl("^Flat", df.out.lst.fg.meanvar$gene.type), ], df.out.lst.bg.meanvar)

barwidth <- 0.8
limits <- aes(ymax = mean.frac + sd.frac, ymin=mean.frac - sd.frac)
ggplot(df.out.lst.meanvar, aes(x = gene.type, y = mean.frac)) + geom_bar(stat = "identity", width = barwidth) + theme_bw() + 
  geom_errorbar(limits, width = barwidth / 2) + theme(aspect.ratio=1.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  
  xlab("") + ylab("Fraction of Genes with Liver-specific DHS peaks")

ggplot(df.out.lst.bg, aes(x = frac.n.spec.by.gene)) + geom_histogram(bins = 25) + 
  theme_bw() + geom_vline(xintercept = c(0.318, 0.225))
