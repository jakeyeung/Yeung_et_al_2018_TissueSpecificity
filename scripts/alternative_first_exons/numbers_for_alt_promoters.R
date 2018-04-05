# 2016-12-30
# Jake Yeung
# # Use Hellinger distance to measure distance between promoter usages
# Copied from alternative_promoter.hellinger.R in liver_kidney_WTKO directory
# Add numbers

rm(list=ls())

setwd("/home/yeung/projects/tissue-specificity")
library(ggrepel)
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

dataset <- "liverWTKO"  # or hogenesch
dataset <- "hogenesch"  # or liverWTKO

distmeth <- "hellinger"
distmeth <- "euclidean"

if (distmeth == "euclidean"){
  jxlab <- "Euclidean Distance"
} else if (distmeth == "hellinger"){
  jxlab <- "Hellinger Distance"
} else {
  stop("distmeth must be Euclidean or Hellinger")
}

do.filter <- FALSE

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
  
  if (do.filter){
    tiss.filt <- c("Liver", "Kidney")
    tpm.afe.avg <- subset(tpm.afe.avg, tissue %in% tiss.filt)
    genes.filt <- as.character(tpm.afe.avg$gene)
    fits.long.filt <- subset(fits.long.filt, gene %in% genes.filt)
  }
  
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
# tpm.afe.dist <- subset(tpm.afe.avg, nprom > 1 & !is.na(amp)) %>%
# tpm.afe.dist <- subset(tpm.afe.avg, nprom > 1) %>%
  group_by(gene) %>%
  do(CalculateDistance(., dist.method = distmeth, jvar = "tpm_norm.avg", transcript_id = "transcript", return.as.df = TRUE))
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

ggplot(tpm.afe.dist, aes(x = proms.dist, y = seq(length(proms.dist)))) + geom_point()


# assigin prom dist to fits long filt
proms.dist.tbl <- hash(as.character(tpm.afe.dist$gene), tpm.afe.dist$proms.dist)

fits.long.filt$proms.dist <- sapply(as.character(fits.long.filt$gene), function(g){
  pdist <- proms.dist.tbl[[g]]
  if (is.null(pdist)){
    pdist <- NA
  } 
  return(pdist)
})

amp.thres <- 2.5 / 2
dist.thres <- 0.55 * sqrt(2)
dist.thres2 <- 0.3 * sqrt(2)
fits.long.filt$label <- mapply(function(amp.avg, prom.dist, gene){
  if (any(is.na(c(amp.avg, prom.dist)))) return(NA)
  if ((amp.avg < amp.thres & prom.dist < dist.thres) | (prom.dist < dist.thres2)){
    return(NA)
  } else {
    return(gene)
  }
}, fits.long.filt$amp.avg, fits.long.filt$proms.dist, as.character(fits.long.filt$gene))

# label Insig2
fits.long.filt$label[[which(fits.long.filt$gene == "Insig2")]] <- "Insig2"

m <- ggplot(fits.long.filt, aes(x = proms.dist / sqrt(2), y = amp.avg* 2, label = label)) + 
  geom_point(alpha = 0.2) + 
  geom_text_repel(size = 5) + 
  # geom_text() + 
  theme_bw(24) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ylab("Average Log2 FC") + xlab(jxlab)
# ylim(0, 5)
print(m)


# Set arbitrary cutoff ----------------------------------------------------

min.euclid.dist <- 0.6

fits.long.filt$has.alt.prom <- sapply(fits.long.filt$proms.dist / sqrt(2), function(d) d >= 0.6)

(n.alt.prom.rhyth <- table(subset(fits.long.filt, select = has.alt.prom)))

(frac.n.alt.prom.rhyth <- n.alt.prom.rhyth[[2]] / sum(n.alt.prom.rhyth))

# How many possible alternative TSS's? ------------------------------------

n.proms <- tpm.afe.avg %>%
  group_by(gene) %>%
  summarise(nprom = unique(nprom))

genes.flat <- as.character(subset(fits.long.filt, model == "")$gene)
genes.ts <- as.character(subset(fits.long.filt, n.rhyth >= 1 & n.rhyth < 8)$gene)

genes.flatts <- hash(c(genes.flat, genes.ts), c(rep("flat", length(genes.flat)), rep("ts", length(genes.ts))))

n.proms$multiprom <- sapply(n.proms$nprom, function(n) n >= 2)

n.proms$rhyth <- sapply(as.character(n.proms$gene), function(g){
  r <- genes.flatts[[g]]
  if (is.null(r)){
    return(NA)
  } else {
    return(r)
  }
})

# make table 
(tbl.rhyth <- table(subset(n.proms, !is.na(rhyth), select = c(rhyth, multiprom))))

# genome wide
(tbl.gw <- table(subset(n.proms, select = c(multiprom))))
n.multiprom <- tbl.gw[[2]] / sum(tbl.gw)  # 0.55


# What fraction of genes are significantly correlated? --------------------


