# 2015-11-19
# Jake Yeung
# Quantify significance of sitecounts.

library(ggplot2)
library(reshape2)
library(dplyr)
library(hash)

dist.ref <- 500  # 500 left and right of promoter is reference
# sitecounts.path <- "/home/yeung/projects/tissue-specificity/data/sitecounts/motevo/sitecount_matrix_geneids"
# dist <- 500
sitecounts.path <- "/home/yeung/projects/tissue-specificity/data/sitecounts/motevo/encode_50000_dist_sum_multigene/sitecounts.50000.multigene.mat"
dist <- 50000  # needs rescaling
# sitecounts.path <- "/home/yeung/projects/tissue-specificity/data/sitecounts/motevo/encode_1000_dist_sum_multigene/sitecounts.1000.multigene.mat"
# dist <- 1000

# Function ----------------------------------------------------------------



# Load --------------------------------------------------------------------

N <- read.table(sitecounts.path, header = TRUE)
load("Robjs/dat.long.fixed_rik_genes.Robj")
load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.20.phase_sd_maxdiff_avg.Robj")
load("Robjs/tpm.gauss.bic_models.Robj")

# ggplot(tpm.gauss, aes(x = center.dists, y = intrascore2, label = gene_name)) + geom_point(alpha = 0.1) + scale_y_log10() + geom_text()
# PromoterSpacePlots.nostics(subset(tpm.gauss, gene_name == "Insig2")$sigs[[1]], "Insig2", draw.ellipse = T)

source("scripts/functions/FisherTestSitecounts.R")
source("scripts/functions/PlotGeneAcrossTissues.R")

# Make long ---------------------------------------------------------------

colnames(N)[1] <- "gene"

N$gene.uniq <- make.names(N$gene, unique = TRUE)

N.long <- melt(N, id.vars = c("gene", "gene.uniq"), variable.name = "motif", value.name = "sitecount")

jgene <- as.character(sample(dat.long$gene, 1))
motifs <- c("RORA.p2", "HNF4A_NR2F1.2.p2", "ONECUT1.2.p2")

N.sub <- subset(N.long, gene == jgene)
for (jmotif in motifs){
  print(subset(N.sub, motif == jmotif))
}


# Enrichment of Liver-specific genes --------------------------------------

key <- as.character(fits.best$gene)
val <- unlist(as.character(fits.best$model))

rhyth.hash <- hash(key, val)

N.long$model <- sapply(N.long$gene, function(g) rhyth.hash[[as.character(g)]])

# do my test
jmotif <- "RORA.p2"
jmotif <- "MEF2.A.B.C.D..p2"
jmotif <- "HNF4A_NR2F1.2.p2"
jmotif <- "ONECUT1.2.p2"
cutoff <- 0.5
N.sub <- subset(N.long, (model == "Liver" | model == "") & motif == jmotif)

N.sub$has.motif <- sapply(N.sub$sitecount, function(s){
  if (s > cutoff){
    return(TRUE)
  } else {
    return(FALSE)
  }
})
N.table <- table(N.sub$has.motif, unlist(N.sub$model))
print(N.table)
print(fisher.test(N.table))



# Enrichment of pairs of motifs -------------------------------------------

scale.factor <- dist / dist.ref

N.sub.base <- subset(N.long, (model == "Liver" | model == ""))
N.sub.base$sitecount.norm <- N.sub.base$sitecount / scale.factor

motifs.all <- sort(unique(as.character(N.long$motif)))
motifs <- c("RORA.p2", sample(motifs.all, 1))
motifs <- c("RORA.p2")
# motifs <- c("MEF2.A.B.C.D..p2")
# motifs <- c("RORA.p2", "HNF4A_NR2F1.2.p2")
# motifs <- c("RORA.p2", "ONECUT1.2.p2")
cutoff <- 0.05

N.sub <- subset(N.sub.base, motif %in% motifs)

N.sub <- N.sub %>%
  group_by(gene.uniq, gene) %>%
  summarise(min.sitecount = min(sitecount.norm))

N.sub$model <- sapply(N.sub$gene, function(g) rhyth.hash[[as.character(g)]])

N.sub$has.motif <- sapply(N.sub$min.sitecount, function(s){
  if (s > cutoff){
    return(TRUE)
  } else {
    return(FALSE)
  }
})
N.table <- table(N.sub$has.motif, unlist(N.sub$model))
print(motifs)
print(N.table)
print(fisher.test(N.table))


# Enrichment for all TFs --------------------------------------------------

N.sub.base <- subset(N.long, (model == "Adr" | model == ""))

start <- Sys.time()
cutoffs <- seq(from = 1, to = 5, by = 0.5)
# cutoffs <- seq(from = 0.5, to = 2.5, by = 0.5)
N.ftest.all <- data.frame()
for (cutoff in cutoffs){
  print(cutoff)
  N.ftest <- N.sub.base %>%
    group_by(motif) %>%
    do(FisherTestSitecounts(., cutoff))
  N.ftest$cutoff <- cutoff
  N.ftest.all <- rbind(N.ftest.all, N.ftest)
}
print(Sys.time() - start)

N.ftest.sum <- N.ftest.all %>%
  group_by(motif) %>%
  summarise(odds.ratio = mean(odds.ratio), p.value = mean(p.value))

ggplot(N.ftest.sum, aes(y = -log10(p.value), x = odds.ratio, label = motif)) + geom_point() + geom_text()

FisherTestSitecounts(subset(N.sub.base, motif == "RORA.p2"), cutoff = 5)

# Enrichment for all TFs: Liver versus [t]issue [w]ide ------------------------------


models.tw <- sort(as.character(unique(subset(fits.best, n.rhyth >= 8)$model)))
# models.tw <- sort(as.character(unique(subset(fits.best, model == "BFAT")$model)))
# models.tw <- sort(as.character(unique(subset(fits.best, model != "Liver" & model != "")$model)))
N.sub.base.livertw <- subset(N.long, model %in% c(models.tw, "Liver"))

length(unique(subset(N.sub.base.livertw, model %in% models.tw)$gene))

N.sub.base.livertw$model <- sapply(N.sub.base.livertw$model, function(m){
  if (m != "Liver"){
    return("TissueWide")
  } else {
    return("Liver")
  }
})

start <- Sys.time()
cutoffs <- seq(from = 1, to = 5, by = 0.5)
# cutoffs <- seq(from = 1, to = 1.5, by = 0.1)
N.ftest.ltw.all <- data.frame()
for (cutoff in cutoffs){
  print(cutoff)
  N.ltw.ftest <- N.sub.base.livertw %>%
    group_by(motif) %>%
    do(FisherTestSitecounts(., cutoff))
  N.ltw.ftest$cutoff <- cutoff
  N.ftest.ltw.all <- rbind(N.ftest.ltw.all, N.ltw.ftest)
}
print(Sys.time() - start)

N.ftest.ltw.sum <- N.ftest.ltw.all %>%
  group_by(motif) %>%
  summarise(odds.ratio = mean(odds.ratio), p.value = mean(p.value))

ggplot(N.ftest.ltw.sum, aes(x = -log10(p.value), y = odds.ratio, label = motif)) + geom_point() + geom_text()


jcutoff <- 2.5
FisherTestSitecounts(subset(N.sub.base.livertw, motif == "ATF5_CREB3.p2"), cutoff=jcutoff, show.table = TRUE)
FisherTestSitecounts(subset(N.sub.base.livertw, motif == "NFKB1_REL_RELA.p2"), cutoff=jcutoff, show.table = TRUE)
FisherTestSitecounts(subset(N.sub.base.livertw, motif == "bHLH_family.p2"), cutoff=jcutoff, show.table = TRUE)
FisherTestSitecounts(subset(N.sub.base.livertw, motif == "RORA.p2"), cutoff=jcutoff, show.table = TRUE)
FisherTestSitecounts(subset(N.sub.base.livertw, motif == "HNF1A.p2"), cutoff=jcutoff, show.table = TRUE)
FisherTestSitecounts(subset(N.sub.base.livertw, motif == "HNF4A_NR2F1.2.p2"), cutoff=jcutoff, show.table = TRUE)
FisherTestSitecounts(subset(N.sub.base.livertw, motif == "MEF2.A.B.C.D..p2"), cutoff=jcutoff, show.table = TRUE)
FisherTestSitecounts(subset(N.sub.base.livertw, motif == "NR6A1.p2"), cutoff=jcutoff, show.table = TRUE)
FisherTestSitecounts(subset(N.sub.base.livertw, motif == "ZBTB16.p2"), cutoff=jcutoff, show.table = TRUE)
FisherTestSitecounts(subset(N.sub.base.livertw, motif == "ONECUT1.2.p2"), cutoff=jcutoff, show.table = TRUE)
FisherTestSitecounts(subset(N.sub.base.livertw, motif == "PITX1..3.p2"), cutoff=jcutoff, show.table = TRUE)

PlotGeneAcrossTissues(subset(dat.long, gene == "Egr1"))

# Compare with promoters --------------------------------------------------

load("Robjs/N.long.promoters_500.Robj")

# models.tw <- sort(as.character(unique(subset(fits.best, n.rhyth >= 8)$model)))
models.tw <- ""

jmodel <- "Mus"
jtiss <- c("Mus")

fits.adrbfataorta <- subset(fits.best, n.rhyth == 3)
fits.adrbfataorta <- fits.adrbfataorta[grep("Adr.*Aorta.*BFAT", fits.adrbfataorta$model), ]
length(unique(fits.adrbfataorta$gene))
jmodel <- unique(as.character(fits.adrbfataorta$model))
jtiss <- jmodel

fits.bfataortamus <- subset(fits.best, n.rhyth == 3)
fits.bfataortamus <- fits.bfataortamus[grep("Aorta.*BFAT.*Mus", fits.bfataortamus$model), ]
jmodel <- unique(as.character(fits.adrbfataorta$model))
jtiss <- jmodel

fits.bfataorta <- subset(fits.best, n.rhyth == 2)
fits.bfataorta <- fits.bfataorta[grep("Aorta.*BFAT|BFAT.*Aorta", fits.bfataorta$model), ]
jmodel <- unique(as.character(fits.bfataorta$model))
# jmodel <- c(jmodel, "BFAT")  # include tissue-specific module BFAT if you want
jtiss <- jmodel

fits.livkid <- subset(fits.best, n.rhyth == 2)
fits.livkid <- fits.livkid[grep("Liver.*Kidney|Kidney.*Liver", fits.livkid$model), ]
jmodel <- unique(as.character(fits.livkid$model))
jtiss <- jmodel

# the Aorta,BFAT antiphasic module
fits.bfataorta <- subset(fits.best, n.rhyth > 1 & n.rhyth < 11)
fits.bfataorta <- fits.bfataorta[grep("(;|^)Aorta.*;BFAT(;|$)", fits.bfataorta$model), ]
jmodel <- unique(as.character(fits.bfataorta$model))
jtiss <- jmodel

# Tissue-wide
fits.tw <- subset(fits.best, n.rhyth >= 8)
jmodel <- unique(as.character(fits.tw$model))
jtiss <- jmodel

N.sub.base.livertw <- subset(N.long, model %in% c(models.tw, jmodel))

length(unique(subset(N.sub.base.livertw, model %in% models.tw)$gene))
# length(fits.livkid <- fits.livkid[grep("Liver.*Kidney|Kidney.*Liver", fits.livkid$model), ]$gene)

N.sub.base.livertw$model <- sapply(N.sub.base.livertw$model, function(m){
  if (!m %in% jtiss){
    return("Flat")
  } else {
    return("Rhyth")
  }
})

length(unique(subset(N.sub.base.livertw, model == "Rhyth")$gene))

start <- Sys.time()
cutoffs <- seq(from = 0.4, to = 0.8, by = 0.1)
N.ftest.ltw.all <- data.frame()
for (cutoff in cutoffs){
  print(cutoff)
  N.ltw.ftest <- N.sub.base.livertw %>%
    group_by(motif) %>%
    do(FisherTestSitecounts(., cutoff))
  N.ltw.ftest$cutoff <- cutoff
  N.ftest.ltw.all <- rbind(N.ftest.ltw.all, N.ltw.ftest)
}
print(Sys.time() - start)

N.ftest.ltw.sum <- N.ftest.ltw.all %>%
  group_by(motif) %>%
  summarise(odds.ratio = mean(odds.ratio), p.value = mean(p.value))

ggplot(N.ftest.ltw.sum, aes(y = -log10(p.value), x = odds.ratio, label = motif)) + geom_point() + geom_text()

jcutoff <- 0.6
FisherTestSitecounts(subset(N.sub.base.livertw, motif == "HNF1A.p2"), cutoff=jcutoff, show.table = TRUE)
FisherTestSitecounts(subset(N.sub.base.livertw, motif == "MEF2.A.B.C.D..p2"), cutoff=jcutoff, show.table = TRUE)
FisherTestSitecounts(subset(N.sub.base.livertw, motif == "NFIL3.p2"), cutoff=jcutoff, show.table = TRUE)
FisherTestSitecounts(subset(N.sub.base.livertw, motif == "RORA.p2"), cutoff=jcutoff, show.table = TRUE)
FisherTestSitecounts(subset(N.sub.base.livertw, motif == "ELK1.4_GABP.A.B1..p3"), cutoff=jcutoff, show.table = TRUE)
FisherTestSitecounts(subset(N.sub.base.livertw, motif == "FOX.C1.C2..p2"), cutoff=jcutoff, show.table = TRUE)
FisherTestSitecounts(subset(N.sub.base.livertw, motif == "bHLH_family.p2"), cutoff=jcutoff, show.table = TRUE)


# Compare with non-expressed genes in Liver -------------------------------

load("Robjs/N.long.dhs_50000_multigene.Robj")

dat.mean <- subset(dat.long, experiment == "rnaseq") %>%
  group_by(gene, tissue) %>%
  summarise(exprs = mean(exprs))

exprs.cutoff <- 2
dat.mean$in.liver <- mapply(function(e, tiss) if (e <= exprs.cutoff & tiss == "Liver") return(TRUE) else return(FALSE), dat.mean$exprs, dat.mean$tissue)

genes.notinliver.all <- subset(dat.mean, in.liver == TRUE)$gene
genes.notinliver <- as.character(intersect(genes.notinliver.all, as.character(fits.best$gene)))  # but expressed in other tissues

genes.liver <- as.character(subset(fits.best, model == "Liver")$gene)

genes.all <- c(genes.liver, genes.notinliver)

N.sub.base.livernotexprs <- subset(N.long, gene %in% genes.all)  # some attrition so ngenes dont always match


N.sub.base.livernotexprs$model <- sapply(N.sub.base.livernotexprs$model, function(m){
  if (m != "Liver"){
    return("Flat")
  } else {
    return("Liver")
  }
})

start <- Sys.time()
cutoffs <- seq(from = 1, to = 5, by = 0.5)
# cutoffs <- seq(from = 1, to = 1.5, by = 0.1)
N.ftest.ltw.all <- data.frame()
for (cutoff in cutoffs){
  print(cutoff)
  N.ltw.ftest <- N.sub.base.livernotexprs %>%
    group_by(motif) %>%
    do(FisherTestSitecounts(., cutoff))
  N.ltw.ftest$cutoff <- cutoff
  N.ftest.ltw.all <- rbind(N.ftest.ltw.all, N.ltw.ftest)
}
print(Sys.time() - start)

N.ftest.ltw.sum <- N.ftest.ltw.all %>%
  group_by(motif) %>%
  summarise(odds.ratio = mean(odds.ratio), p.value.minuslog = mean(-log10(p.value)))

ggplot(N.ftest.ltw.sum, aes(x = p.value.minuslog, y = odds.ratio, label = motif)) + geom_point() + geom_text()


jcutoff <- 5
FisherTestSitecounts(subset(N.sub.base.livernotexprs, motif == "ONECUT1.2.p2"), cutoff=jcutoff, show.table = TRUE)
FisherTestSitecounts(subset(N.sub.base.livernotexprs, motif == "HNF1A.p2"), cutoff=jcutoff, show.table = TRUE)
FisherTestSitecounts(subset(N.sub.base.livernotexprs, motif == "HNF4A_NR2F1.2.p2"), cutoff=jcutoff, show.table = TRUE)
FisherTestSitecounts(subset(N.sub.base.livernotexprs, motif == "NR6A1.p2"), cutoff=jcutoff, show.table = TRUE)
FisherTestSitecounts(subset(N.sub.base.livernotexprs, motif == "PITX1..3.p2"), cutoff=jcutoff, show.table = TRUE)

PlotGeneAcrossTissues(subset(dat.long, gene == "Pitx3"))
