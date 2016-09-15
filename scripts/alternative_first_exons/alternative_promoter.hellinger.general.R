# 2016-09-15
# Jake Yeung
# # Use Hellinger distance to measure distance between promoter usages
# Copied from alternative_promoter.hellinger.R in liver_kidney_WTKO directory

rm(list=ls())

setwd("/home/yeung/projects/tissue-specificity")

library(dplyr)
library(ggplot2)
library(hash)
library(reshape2)

source("scripts/functions/BiomartFunctions.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/LiverKidneyFunctions.R")
source("scripts/functions/AlternativeFirstExonsFunctions.R")
source("scripts/functions/NcondsFunctions.R")
source("scripts/functions/ModelStrToModel.R")

eps <- 1  # for log2 transform


# Load --------------------------------------------------------------------

dataset <- "hogenesch"  # or liverWTKO
dataset <- "liverWTKO"  # or hogenesch

if (dataset == "hogenesch"){
  load("Robjs/tpm.afe.avg.binary.Robj", verbose=T)
  load("Robjs/dat.long.fixed_rik_genes.Robj", verbose=T)
  load("Robjs/tpm.merged.Robj", verbose=T); dat.bytranscript <- tpm.merged; rm(tpm.merged)
  load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.20.phase_sd_maxdiff_avg.Robj", v=T); fits.long.filt <- fits.best; rm(fits.best)
  tpm.afe.avg <- dplyr::rename(tpm.afe.avg, gene = gene_name)
  tpm.afe.avg <- dplyr::rename(tpm.afe.avg, transcript = transcript_id)
  tpm.afe.avg <- dplyr::rename(tpm.afe.avg, tpm = tpm.avg)
  
  tpm.counts <- subset(tpm.afe.avg, tissue == "Liver") %>%
    group_by(tissue, gene) %>%
    summarise(counts = length(transcript))
  tpm.counts <- tpm.counts[order(tpm.counts$counts, decreasing = T), ]
  
  counts.dic <- hash(as.character(tpm.counts$gene), tpm.counts$counts)
  
  # tiss.filt <- c("Liver", "Kidney")
  # tpm.afe.avg <- subset(tpm.afe.avg, tissue %in% tiss.filt)
  
} else if (dataset == "liverWTKO"){
  load("Robjs/liver_kidney_atger_nestle/tpm.afe.avg.binary.Robj", v=T)
  load("Robjs/liver_kidney_atger_nestle/dat.long.liverkidneyWTKO.bugfixed.Robj", v=T)
  load("Robjs/liver_kidney_atger_nestle/dat.bytranscript.bugfixed.Robj", v=T)
  load("Robjs/liver_kidney_atger_nestle/fits.long.multimethod.filtbest.staggeredtimepts.bugfixed.annotated.Robj", v=T)
  dat.bytranscript <- StaggeredTimepointsLivKid(dat.bytranscript)
  dat.bytranscript <- CollapseTissueGeno(dat.bytranscript)  # match fits.long.filt
  dat.long <- StaggeredTimepointsLivKid(dat.long)
  fits.long.filt <- subset(fits.long.filt, method == "g=1001")
} else {
  stop("Dataset mustt be hogenesch or liverWTKO")
}
print(head(tpm.afe.avg))

# Load data by transcript


# do only on subset of models??
# jgenes <- unique(as.character(subset(fits.long.filt, model != "")$gene))

# Calculate Hellinger distance  -------------------------------------------

jtest <- subset(tpm.afe.avg, nprom > 1 & gene %in% c("Insig2", "Upp2", "Slc45a3", "Ddc"))
jtest <- subset(tpm.afe.avg, nprom > 1 & gene %in% c("Insig2"))
jtest <- subset(tpm.afe.avg, nprom > 1 & gene %in% c("0610007P14Rik"))
jtest <- subset(tpm.afe.avg, nprom > 1 & gene %in% c("6430550D23Rik"))
jtest <- subset(tpm.afe.avg, nprom > 1 & gene %in% c("Upp2"))
jtest <- subset(tpm.afe.avg, nprom > 1 & gene %in% c("Ddc"))
jtest$tissue <- ifelse(jtest$amp, yes = "rhyth.tiss", no = "flat.tiss")

start <- Sys.time()

fits.long.filt$n.tiss <- sapply(as.character(fits.long.filt$model), function(m) return(length(ModelToTissue(m))))
max.tiss <- max(fits.long.filt$n.tiss)
genes.to.check <- as.character(subset(fits.long.filt, n.tiss > 0 & n.tiss < max.tiss)$gene)

tpm.afe.avg$tissue <- ifelse(tpm.afe.avg$amp, yes = "rhyth.tiss", no = "flat.tiss")
print(subset(tpm.afe.avg, gene == "Insig2"))
tpm.afe.avg <- tpm.afe.avg %>%
  group_by(gene, transcript, tissue, nprom, amp) %>%
  summarise(tpm = sum(tpm)) %>%
  group_by(gene, tissue) %>%
  mutate(tpm_norm.avg = CalculateFractionIsoformUsage(tpm, pseudocount = 1))
print(subset(tpm.afe.avg, gene == "Insig2"))


tpm.afe.dist <- subset(tpm.afe.avg, nprom > 1 & !is.na(amp) & gene %in% genes.to.check) %>%
  group_by(gene) %>%
  do(CalculateDistance(., dist.method = "hellinger", jvar = "tpm_norm.avg", transcript_id = "transcript", return.as.df = TRUE))
print(Sys.time() - start)

gene.models <- hash(as.character(fits.long.filt$gene), as.character(fits.long.filt$model))
tpm.afe.dist$model <- sapply(as.character(tpm.afe.dist$gene), function(g) gene.models[[g]])
# take top N models
top.n <- 15
tpm.model.count <- tpm.afe.dist %>%
  group_by(model) %>%
  summarise(count = length(gene)) %>%
  arrange(desc(count))
top.models <- tpm.model.count$model[1:top.n]
jsub <- subset(tpm.afe.dist, model %in% top.models)
ggplot(jsub, aes(x = model, y = proms.dist)) + geom_boxplot()
  
# 
# # Compare rhythmic genes with flat genes ----------------------------------
# 
# gene.models <- hash(as.character(fits.long.filt$gene), as.character(fits.long.filt$model))

# tpm.afe.dist$model <- sapply(as.character(tpm.afe.dist$gene), function(g) gene.models[[g]])

# 
# # show only models with >150 genes
# models.count <- fits.long.filt %>%
#   group_by(model) %>%
#   summarise(count = length(gene))
# signif.models <- as.character(subset(models.count, count > 150)$model)
# 
# # filter genes similarly expressed
# dat.diff <- subset(dat.long) %>%
#   group_by(gene, tissue) %>%
#   summarise(exprs.mean = mean(exprs)) %>%
#   group_by(gene) %>%
#   summarise(exprs.diff = diff(exprs.mean))
# jcutoff <- Inf
# jgenes.similar <- as.character(subset(dat.diff, abs(exprs.diff) <= jcutoff)$gene)
# 
# # signif.models <- c("", fg.model)
# signif.models <- c("", "Liver_SV129", "Liver_SV129,Liver_BmalKO", "Liver_SV129,Kidney_SV129", "Kidney_SV129", "Kidney_SV129,Kidney_BmalKO")
# # signif.models <- c("", "Liver_SV129")
# # signif.models <- c("", "Liver_SV129,Liver_BmalKO")
# jsub <- subset(tpm.afe.dist, model %in% signif.models & gene %in% jgenes.similar)
# ggplot(jsub, aes(x = model, y = proms.dist)) + geom_boxplot()
# 
# fg.model <- "Liver_SV129"
# fg.model <- "Kidney_SV129,Kidney_BmalKO"
# bg.model <- ""
# jsub2.flat <- subset(jsub, model %in% c(bg.model))
# jsub2.fg <- subset(jsub, model %in% c(fg.model))
# ks.test(jsub2.flat$proms.dist, jsub2.fg$proms.dist)
# t.test(jsub2.flat$proms.dist, jsub2.fg$proms.dist)
# 
# 
# 
# # Investigate why Flat genes have high alt promoter usage -----------------
# 
# jorder <- jsub[order(jsub$proms.dist, decreasing = TRUE), ]
# jgene <- "Rpl13a"
# jgene <- "Uba52"
# jgene <- "Slc45a3"
# jgene <- "Birc6"
# jgene <- "Uba52"
# jgene <- "Gm14327"
# 
# PlotGeneTissuesWTKO(subset(dat.long, gene == jgene))
# proms.sub <- GetPromoterUsage(subset(tpm.afe.avg, gene == jgene), do.svd = FALSE, transcript_id = "transcript", get.prom.only = TRUE)
# # proms.sub.bysamp <- GetPromoterUsage(subset(tpm.afe.bysamp, gene == jgene & time == 22), do.svd = FALSE, transcript_id = "transcript", get.prom.only = TRUE)
# 
# print(PlotTpmAcrossTissues(subset(dat.bytranscript, gene == jgene & geno == "SV129"), jtitle = jgene, log2.transform = TRUE, transcript_id = "transcript") + facet_wrap(~tissue))
# print(jgene)
# 
# proms.sub <- GetPromoterUsage(subset(tpm.afe.avg, gene == jgene), do.svd = FALSE, transcript_id = "transcript", get.prom.only = TRUE)
# print(dist(proms.sub))
# proms.sub.bysamp <- GetPromoterUsage(subset(tpm.afe.bysamp, gene == jgene), do.svd = FALSE, transcript_id = "transcript", get.prom.only = TRUE)
# print(dist(proms.sub.bysamp))
# 
# dist.by.time <- subset(tpm.afe.bysamp, gene == jgene) %>%
#   group_by(time) %>%
#   do(jdist = as.numeric(dist(GetPromoterUsage(., do.svd=FALSE, transcript_id="transcript", get.prom.only=TRUE))))
# dist.by.time$dist <- sapply(dist.by.time$jdist, function(s) s[[1]])
# ggplot(dist.by.time, aes(x = time, y = dist)) + geom_bar(stat = "identity")
# 
# jmodel <- as.character(subset(fits.long.filt, gene == jgene)$model[[1]])
# 
# tx <- GetPromoterUsage(subset(tpm.afe.bysamp, gene == jgene & time == 2), transcript_id = "transcript", get.entropy=FALSE, return.transcripts=TRUE)[[1]]$transcript
# print(PlotGeneTissuesWTKO(subset(dat.long, gene == jgene)) + ggtitle(paste(jgene, jmodel)))
# PromoterSpacePlots(tpm.afe.bysamp, jgene = jgene, transcript_id = "transcript", use.weights = FALSE)
# print(PlotTpmAcrossTissues(subset(dat.bytranscript, gene == jgene & transcript %in% tx), jtitle = jgene, log2.transform = TRUE, transcript_id = "transcript"))
# print(PlotTpmAcrossTissues(subset(dat.bytranscript, gene == jgene), jtitle = jgene, log2.transform = TRUE, transcript_id = "transcript"))
# 
