# 2016-07-08
# Jake Yeung
# Do alt promoter per sample


rm(list=ls())

library(dplyr)
library(ggplot2)
library(hash)
library(mvtnorm)
library(reshape2)
library(wordcloud)

source("scripts/functions/BiomartFunctions.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/LiverKidneyFunctions.R")
source("scripts/functions/AlternativeFirstExonsFunctions.R")
source("scripts/functions/NcondsFunctions.R")

eps <- 1  # for log2 transform

DifferentialPromoterUsageTwoTissues <- function(dat.sub, show.plot=FALSE, jvar = "tpm.norm", component = 1){
  # dat.sub <- subset(dat.bytranscript, gene == jgene & geno == jgeno)
  jsub.mat <- dcast(dat.sub, formula = transcript ~ tissue + time + geno, value.var = jvar)
  rownames(jsub.mat) <- jsub.mat$transcript; jsub.mat$transcript <- NULL
  
  # print(as.character(dat.sub$gene[[1]]))
  # print(jsub.mat)

  if (any(is.nan(as.matrix(jsub.mat)))){
    out.df <- data.frame(NULL)
    return(out.df)
  }
  s <- svd(jsub.mat)
  
  rownames(s$v) <- colnames(jsub.mat)
  
  # do ttest
  if (show.plot){
    tissues <- sapply(colnames(jsub.mat), function(cname) strsplit(cname, "_")[[1]][[1]], USE.NAMES = FALSE)
    jtiss <- unique(tissues)
    jgene <- dat.sub$gene[[1]]
    boxplot(s$v[, 1][grepl(jtiss[[1]], rownames(s$v))], s$v[, 1][grepl(jtiss[[2]], rownames(s$v))], names = jtiss, main = jgene)
  }
  
  tissues <- sapply(colnames(jsub.mat), function(cname) strsplit(cname, "_")[[1]][[1]], USE.NAMES = FALSE)
  eigensamp <- s$v[, component] 
  # eigenpromoter <- s$u[, component] * s$d[[component]]
  
  out = tryCatch({
    t.test(eigensamp ~ tissues)
  }, error = function(e) {
    list()
  })
  
  # out <- t.test(eigensamp ~ tissues)
  out.df <- data.frame(mean1 = out$estimate[[1]], mean2 = out$estimate[[2]], pval = out$p.value)
  return(out.df)
}

GetTranscriptLoadings <- function(dat.sub, component = 1, jvar = "tpm.norm"){
  jsub.mat <- dcast(dat.sub, formula = transcript ~ tissue + time + geno, value.var = jvar)
  rownames(jsub.mat) <- jsub.mat$transcript; jsub.mat$transcript <- NULL
  s <- svd(jsub.mat)
  
  rownames(s$u) <- rownames(jsub.mat)
  return(s$u[, 1])
}

# Load --------------------------------------------------------------------

# load("Robjs/liver_kidney_atger_nestle/dat.bytranscript.bugfixed.Robj", v=T)
load("/home/yeung/projects/tissue-specificity/Robjs/liver_kidney_atger_nestle/dat.bytranscript.bugfixed.annotated.Robj", v=T)
load("Robjs/liver_kidney_atger_nestle/dat.long.liverkidneyWTKO.bugfixed.Robj", v=T)
load("Robjs/liver_kidney_atger_nestle/fits.long.multimethod.filtbest.staggeredtimepts.bugfixed.Robj")

dat.bytranscript <- StaggeredTimepointsLivKid(dat.bytranscript)
dat.bytranscript <- CollapseTissueGeno(dat.bytranscript)  # match fits.long.filt

dat.orig <- dat.long
dat.long <- StaggeredTimepointsLivKid(dat.long)

# Calculate fractional isoform usage --------------------------------------

dat.bytranscript <- dat.bytranscript %>%
  group_by(gene, tissue, time, geno) %>%
  mutate(tpm.norm = tpm / sum(tpm))

# label npromoters
tpm.counts <- subset(dat.bytranscript, tissue == "Liver_SV129" & time == 2) %>%
  group_by(tissue, gene) %>%
  summarise(counts = length(transcript))
tpm.counts <- tpm.counts[order(tpm.counts$counts, decreasing = T), ]
counts.dic <- hash(as.character(tpm.counts$gene), tpm.counts$counts)
dat.bytranscript$nprom <- sapply(as.character(dat.bytranscript$gene), function(jgene) counts.dic[[jgene]])

# get mean exprs
dat.mean <- subset(dat.long, geno == "SV129" & gene != "") %>%
  group_by(gene) %>%
  summarise(exprs.mean = mean(exprs))

exprs.hash <- hash(as.character(dat.mean$gene), dat.mean$exprs.mean)
dat.bytranscript$exprs.mean <- sapply(as.character(dat.bytranscript$gene), function(k){
  exprs.mean <- exprs.hash[[k]]
  if (is.null(exprs.mean)){
    exprs.mean <- 0
  }
  return(exprs.mean)
})

jgeno <- "SV129"
nprom.cutoff <- 2  # or more
exprs.cutoff <- 2  # log2
dat.altpromtest <- subset(dat.bytranscript, geno == jgeno & nprom >= nprom.cutoff & exprs.mean > 2) %>%
  group_by(gene) %>% 
  do(DifferentialPromoterUsageTwoTissues(., show.plot=FALSE, jvar = "tpm.norm"))


# Sanity checks -----------------------------------------------------------

jmods <- c("Liver_SV129", "Kidney_SV129")
jgenes <- unique(as.character(subset(fits.long.filt, model %in% jmods)$gene))

dat.altfilt <- subset(dat.altpromtest, gene %in% jgenes)

dat.altfilt <- dat.altfilt[order(dat.altfilt$pval), ]

genes <- as.character(head(dat.altfilt, n = 50)$gene)

pdf("plots/alternative_exon_usage/liver_kidney.atger_nestle.ttest.pdf")
for (jgene in genes){
  print(jgene)
  tx <- GetTranscriptLoadings(subset(dat.bytranscript, gene == jgene & geno == "SV129"))
  tx <- names(tx[which.max(abs(tx))])
  DifferentialPromoterUsageTwoTissues(subset(dat.bytranscript, gene == jgene & geno == "SV129"), show.plot = TRUE)
  print(PlotGeneTissuesWTKO(subset(dat.long, gene == jgene)) + ggtitle(jgene))
  print(PlotTpmAcrossTissues(subset(dat.bytranscript, gene == jgene & geno == "SV129" & transcript %in% tx), jtitle = jgene, log2.transform = TRUE, transcript_id = "transcript"))
  print(PlotTpmAcrossTissues(subset(dat.bytranscript, gene == jgene & geno == "SV129"), jtitle = jgene, log2.transform = TRUE, transcript_id = "transcript"))
}
dev.off()

# 
# jgene <- "Psen2"
# jgene <- "Ddc"
# jgene <- "Slc45a3"
# jgene <- "Insig2"
# jgene <- "1110002E22Rik"
# 
# out.df <- DifferentialPromoterUsageTwoTissues(subset(dat.bytranscript, gene == jgene & geno == "SV129"), show.plot=TRUE, jvar = "tpm.norm")
# 
# jgeno <- "SV129"
# jtiss <- c("Kidney", "Liver")
# jsub <- subset(dat.bytranscript, gene == jgene & geno == jgeno)
# jsub.mat <- dcast(jsub, formula = transcript ~ tissue + time + geno)
# rownames(jsub.mat) <- jsub.mat$transcript; jsub.mat$transcript <- NULL
# s <- svd(jsub.mat)
# 
# rownames(s$v) <- colnames(jsub.mat)
# rownames(s$u) <- rownames(jsub.mat)
# 
# # do ttest
# tissues <- sapply(colnames(jsub.mat), function(cname) strsplit(cname, "_")[[1]][[1]], USE.NAMES = FALSE)
# out <- t.test(s$v[, 1] ~ tissues)
# 
# boxplot(s$v[, 1][grepl(jtiss[[1]], rownames(s$v))], s$v[, 1][grepl(jtiss[[2]], rownames(s$v))], names = jtiss)
# 
# # plot(s$v[, 1] * s$d[[1]], s$v[, 2] * s$d[[2]])
# # text(s$v[, 1] * s$d[[1]], s$v[, 2] * s$d[[2]], colnames(jsub.mat))
# # GetPromoterUsage(jsub)
