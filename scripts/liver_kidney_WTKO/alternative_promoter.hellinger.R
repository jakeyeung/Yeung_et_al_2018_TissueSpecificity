# 2016-08-29
# Jake Yeung
# # Use Hellinger distance to measure distance between promoter usages

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

eps <- 1  # for log2 transform


# Load --------------------------------------------------------------------

load("Robjs/liver_kidney_atger_nestle/dat.bytranscript.bugfixed.Robj", v=T)
load("Robjs/liver_kidney_atger_nestle/dat.long.liverkidneyWTKO.bugfixed.Robj", v=T)
load("Robjs/liver_kidney_atger_nestle/fits.long.multimethod.filtbest.staggeredtimepts.bugfixed.Robj", v=T)
load("Robjs/liver_kidney_atger_nestle/fits.long.multimethod.filtbest.staggeredtimepts.bugfixed.annotated.Robj", v=T)

fits.orig <- fits.long.filt
fits.long.filt <- subset(fits.long.filt, method == "g=1001")

if (!is.null(fits.long.filt$n.params)){
  # Annotate fits
  fits.long.filt$n.params <- sapply(fits.long.filt$model, function(m) return(length(strsplit(as.character(m), ";")[[1]])))
  fits.long.filt$n.rhyth <- sapply(fits.long.filt$model, GetNrhythFromModel)
  fits.long.filt$amp.avg <- mapply(GetAvgAmpFromParams, fits.long.filt$param.list, fits.long.filt$model)
  fits.long.filt$phase.sd <- mapply(GetSdPhaseFromParams, fits.long.filt$param.list, fits.long.filt$model)
  fits.long.filt$phase.maxdiff <- mapply(GetMaxPhaseDiffFromParams, fits.long.filt$param.list, fits.long.filt$model)
  fits.long.filt$phase.avg <- mapply(GetPhaseFromParams, fits.long.filt$param.list, fits.long.filt$model)
}

dat.bytranscript <- StaggeredTimepointsLivKid(dat.bytranscript)
dat.bytranscript <- CollapseTissueGeno(dat.bytranscript)  # match fits.long.filt

dat.long <- StaggeredTimepointsLivKid(dat.long)

# load("Robjs/liver_kidney_atger_nestle/tpm.afe.avg.binary.Robj", v=T)

if (!exists("tpm.afe.avg")){
  tpm.afe.avg <- dat.bytranscript %>%
    mutate(tiss.temp = strsplit(as.character(tissue), "_")[[1]][[1]]) %>%  # Liver_SV129 -> Liver
    group_by(gene, transcript, tiss.temp) %>%
    summarise(tpm = sum(tpm)) %>%
    group_by(gene, tiss.temp) %>%
    mutate(tpm_norm.avg = CalculateFractionIsoformUsage(tpm, pseudocount = 1))
  tpm.afe.avg$tissue <- tpm.afe.avg$tiss.temp
  
  # Count promoters ---------------------------------------------------------
  tpm.counts <- subset(tpm.afe.avg, tissue == "Liver") %>%
    group_by(tissue, gene) %>%
    summarise(counts = length(transcript))
  tpm.counts <- tpm.counts[order(tpm.counts$counts, decreasing = T), ]
  
  counts.dic <- hash(as.character(tpm.counts$gene), tpm.counts$counts)
  tpm.afe.avg$nprom <- sapply(as.character(tpm.afe.avg$gene), function(jgene) counts.dic[[jgene]])
  
  # hash it
  tiss <- c("Liver", "Kidney")
  keys.df <- fits.long.filt %>%
    group_by(gene) %>%
    do(GetGeneModelKeys(., tiss))
  keys <- paste(keys.df$tissue, keys.df$gene, sep = ",")
  vals <- as.numeric(keys.df$is.rhyth)
  is.rhyth.dic <- hash(keys, vals)
  
  # annotate it
  
  # 1 or 0 to signify whether it is rhythmic
  tpm.afe.avg$amp <- mapply(function(jtiss, jgene){
    a <- is.rhyth.dic[[paste(jtiss, jgene, sep = ",")]]
    if (is.null(a)){
      a <- NA
    } 
    return(a)
  }, as.character(tpm.afe.avg$tissue), as.character(tpm.afe.avg$gene))
  outf <- "Robjs/liver_kidney_atger_nestle/tpm.afe.avg.binary.Robj"
  if (!file.exists(outf)){
    save(tpm.afe.avg, file = outf)
  } else {
    print("Skipping saving of Robj")
  }
}

# do only on subset of models??
jgenes <- unique(as.character(subset(fits.long.filt, model != "")$gene))

# Calculate Hellinger distance  -------------------------------------------

jtest <- subset(tpm.afe.avg, nprom > 1 & gene %in% c("Insig2", "Upp2", "Slc45a3", "Ddc"))
start <- Sys.time()
tpm.afe.dist <- subset(tpm.afe.avg, nprom > 1 & !is.na(amp)) %>%
  group_by(gene) %>%
  do(CalculateDistance(., dist.method = "hellinger", jvar = "tpm_norm.avg", transcript_id = "transcript", return.as.df = TRUE))
print(Sys.time() - start)
  

# Compare rhythmic genes with flat genes ----------------------------------

gene.models <- hash(as.character(fits.long.filt$gene), as.character(fits.long.filt$model))

tpm.afe.dist$model <- sapply(as.character(tpm.afe.dist$gene), function(g) gene.models[[g]])

# show only models with >150 genes
models.count <- fits.long.filt %>%
  group_by(model) %>%
  summarise(count = length(gene))
signif.models <- as.character(subset(models.count, count > 150)$model)

# filter genes similarly expressed
dat.diff <- subset(dat.long) %>%
  group_by(gene, tissue) %>%
  summarise(exprs.mean = mean(exprs)) %>%
  group_by(gene) %>%
  summarise(exprs.diff = diff(exprs.mean))
jcutoff <- Inf
jgenes.similar <- as.character(subset(dat.diff, abs(exprs.diff) <= jcutoff)$gene)

# signif.models <- c("", fg.model)
signif.models <- c("", "Liver_SV129", "Liver_SV129,Liver_BmalKO", "Liver_SV129,Kidney_SV129", "Kidney_SV129", "Kidney_SV129,Kidney_BmalKO")
# signif.models <- c("", "Liver_SV129")
# signif.models <- c("", "Liver_SV129,Liver_BmalKO")
jsub <- subset(tpm.afe.dist, model %in% signif.models & gene %in% jgenes.similar)
ggplot(jsub, aes(x = model, y = proms.dist)) + geom_boxplot()

fg.model <- "Liver_SV129"
fg.model <- "Kidney_SV129,Kidney_BmalKO"
bg.model <- ""
jsub2.flat <- subset(jsub, model %in% c(bg.model))
jsub2.fg <- subset(jsub, model %in% c(fg.model))
ks.test(jsub2.flat$proms.dist, jsub2.fg$proms.dist)
t.test(jsub2.flat$proms.dist, jsub2.fg$proms.dist)



# Investigate why Flat genes have high alt promoter usage -----------------

jorder <- jsub[order(jsub$proms.dist, decreasing = TRUE), ]
jgene <- "Rpl13a"
jgene <- "Uba52"
jgene <- "Slc45a3"
jgene <- "Birc6"
jgene <- "Uba52"
jgene <- "Gm14327"

PlotGeneTissuesWTKO(subset(dat.long, gene == jgene))
proms.sub <- GetPromoterUsage(subset(tpm.afe.avg, gene == jgene), do.svd = FALSE, transcript_id = "transcript", get.prom.only = TRUE)
# proms.sub.bysamp <- GetPromoterUsage(subset(tpm.afe.bysamp, gene == jgene & time == 22), do.svd = FALSE, transcript_id = "transcript", get.prom.only = TRUE)

print(PlotTpmAcrossTissues(subset(dat.bytranscript, gene == jgene & geno == "SV129"), jtitle = jgene, log2.transform = TRUE, transcript_id = "transcript") + facet_wrap(~tissue))
print(jgene)

proms.sub <- GetPromoterUsage(subset(tpm.afe.avg, gene == jgene), do.svd = FALSE, transcript_id = "transcript", get.prom.only = TRUE)
print(dist(proms.sub))
proms.sub.bysamp <- GetPromoterUsage(subset(tpm.afe.bysamp, gene == jgene), do.svd = FALSE, transcript_id = "transcript", get.prom.only = TRUE)
print(dist(proms.sub.bysamp))

dist.by.time <- subset(tpm.afe.bysamp, gene == jgene) %>%
  group_by(time) %>%
  do(jdist = as.numeric(dist(GetPromoterUsage(., do.svd=FALSE, transcript_id="transcript", get.prom.only=TRUE))))
dist.by.time$dist <- sapply(dist.by.time$jdist, function(s) s[[1]])
ggplot(dist.by.time, aes(x = time, y = dist)) + geom_bar(stat = "identity")

jmodel <- as.character(subset(fits.long.filt, gene == jgene)$model[[1]])

tx <- GetPromoterUsage(subset(tpm.afe.bysamp, gene == jgene & time == 2), transcript_id = "transcript", get.entropy=FALSE, return.transcripts=TRUE)[[1]]$transcript
print(PlotGeneTissuesWTKO(subset(dat.long, gene == jgene)) + ggtitle(paste(jgene, jmodel)))
PromoterSpacePlots(tpm.afe.bysamp, jgene = jgene, transcript_id = "transcript", use.weights = FALSE)
print(PlotTpmAcrossTissues(subset(dat.bytranscript, gene == jgene & transcript %in% tx), jtitle = jgene, log2.transform = TRUE, transcript_id = "transcript"))
print(PlotTpmAcrossTissues(subset(dat.bytranscript, gene == jgene), jtitle = jgene, log2.transform = TRUE, transcript_id = "transcript"))

