# 2016-07-13
# Jake Yeung
# Do analysis on restricted feeding to check
# Kallisto outliers?

library(dplyr)
library(ggplot2)
library(reshape2)
library(PhaseHSV)
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/LiverKidneyFunctions.R")
source("scripts/functions/BiomartFunctions.R")
source("scripts/functions/PlotFunctions.R")


# Load Kallisto -----------------------------------------------------------

load("Robjs/liver_kidney_atger_nestle/dat.long.liverkidneyWTKO.bugfixed.Robj", v=T)
dat.orig <- dat.long
dat.long <- CollapseTissueGeno(dat.long, keep.tissue.col = TRUE)
dat.long <- StaggeredTimepointsLivKid(dat.long)


# Load Cedric Liver and Kidney WTKO ---------------------------------------

dat.livkid <- LoadLivKid()


# Load Exon-Intron quantification -----------------------------------------

inf <- "/home/yeung/data/tissue_specificity/atger_et_al/GSE73554_WT_RF_Intron_Exon_RFP.txt"
# inf <- "/home/yeung/data/tissue_specificity/atger_et_al/GSE73554_WT_AL_Intron_Exon_RFP.txt"
dat <- read.delim2(inf, header = TRUE, sep = "\t")

keep.cols <- "*Exon*"

dat.filt <- dat[, grepl(keep.cols, colnames(dat)), ]

rownames(dat.filt) <- make.unique(as.character(dat$Gene_Symbol))

times <- as.numeric(sapply(colnames(dat.filt), function(cname) strsplit(cname, "_")[[1]][[4]]))
samps <- sapply(colnames(dat.filt), function(cname) strsplit(cname, "_")[[1]][[5]])

dat.atger.long <- data.frame(gene = dat$Gene_Symbol,
                             tissue = "Liver",
                             time = rep(times, each = nrow(dat.filt)),
                             experiment = "atger",
                             exprs = as.numeric(as.character(unlist(dat.filt))),
                             samp = rep(samps, each = nrow(dat.filt)))

dat.atger.long$samp.numeric <- as.numeric(chartr(old = "CDAB", new = '1234', x = dat.atger.long$samp))

# dat.atger.long <- subset(dat.atger.long, samp %in% c("A", "B"))
dat.atger.long <- subset(dat.atger.long, samp %in% c("C", "D"))
dat.atger.long$time <- mapply(function(time, samp.num) time + 24 * (samp.num - 1), dat.atger.long$time, dat.atger.long$samp.numeric)

# samples C and D are SV129, A and B are C57B6
jgene <- "Per1"
jgene <- "Npas2"
jgene <- "Nfil3"
jgene <- "Nfasc"
jgene <- "Zfp618"
jgene <- "Igkj2"
jgene <- "Arntl"
jgene <- "Ppard"
jgene <- "Snurf"
jgene <- "Ighj1"
jgene <- "Gm8325"
jgene <- "Igkj2"
jgene <- "Gm26602"
jgene <- "Il6st"
jgene <- "Wrnip1"
jgene <- "Arntl"
jgene <- "Gm17181"
jgene <- "Mt1"
jgene <- "Mreg"
jgene <- "Upp2"
jgene <- "Wrnipl1"
jgene <- "Dbp"
jgene <- "Cyp2a4"
jgene <- "Hspa1b"
jgene <- "Pfkfb3"
jgene <- "Mfsd2a"
jgene <- "Gm17181"
jgene <- "Akr1c18"
jgene <- "Sdf2l1"
jgene <- "Gm8325"
PlotGeneTissuesWTKO(subset(dat.orig, gene == jgene)) + theme_bw() + theme(aspect.ratio = 1) + ggtitle(jgene)
m1 <- PlotGeneAcrossTissues(subset(dat.atger.long, gene == jgene)) + theme_bw() + theme(aspect.ratio = 1)
m2 <- PlotGeneAcrossTissues(subset(dat.long, gene == jgene & tissue == "Liver_SV129")) + theme_bw() + theme(aspect.ratio = 1)
multiplot(m1, m2, cols = 2)


  # do pca

# PCA ---------------------------------------------------------------------

dat.type <- "kallisto"  # or kallisto exon livkid orig
if (dat.type == "exon"){
  dat.atger.mean <- dat.atger.long %>%
    group_by(gene, tissue, time, samp) %>%
    summarise(exprs = mean(exprs))
  M <- dcast(dat.atger.mean, formula = gene ~ time, value.var = "exprs")
} else if (dat.type == "kallisto"){
  # filter for protein coding only??
  # jtiss <- "Kidney_SV129"
  jtiss <- "Liver_SV129"
  genes.all <- as.character(unique(dat.long$gene))
  genes.all <- genes.all[which(!is.na(genes.all))]
  genes.all.status <- AnnotatePseudogenes(genes.all, return.original = TRUE)
  genes.hash <- hash(genes.all, genes.all.status)
  dat.long$status <- sapply(as.character(dat.long$gene), function(g){
    status <- genes.hash[[g]]
    if (is.null(status)){
      return(NA)
    } else {
      return(status)
    }
  })
  dat.long.sub <- subset(dat.long, status == "protein_coding")
  M <- dcast(subset(dat.long.sub, tissue == jtiss & !is.na(gene)), formula = gene ~ time, value.var = "exprs")
} else if (dat.type == "livkid"){
  M <- dcast(subset(dat.livkid, tissue == "Liver"), formula = gene ~ time , value.var = "exprs")
  # M <- dcast(subset(dat.livkid, tissue == "Kidney"), formula = gene ~ time , value.var = "exprs")
} else if (dat.type == "orig"){
  # if (!exists(dat.orig.rm)){
  #   dat.orig.rm <- RemoveLowExprsPseudoShortGenes(dat.orig)
  # }
  dat.orig.rm <- dat.orig
  jtiss <- "Kidney"
  M <- dcast(subset(dat.orig.rm, tissue == jtiss & geno == "SV129" & !is.na(gene)), formula = gene ~ time, value.var = "exprs")
}
rownames(M) <- M$gene; M$gene <- NULL

labs <- colnames(M)
phases <- as.numeric(sapply(labs, function(l) strsplit(l, "_")[[1]][[1]]))
phases <- sapply(phases, function(p){
  if (p > 48) return(p - 48)
  if (p > 24) return(p - 24)
  return(p)
})
phases.rad <- phases * 2 * pi / 24
cols <- hsv(PhaseToHsv(phases.rad, min.phase = 0, max.phase = 2 * pi), s = 1, v = 1)

# row center
M.centered <- t(scale(t(M), center = TRUE, scale = FALSE))

dat.pca <- prcomp(M.centered, center = F, scale. = F)

pc1 <- 1
pc2 <- 2
plot(dat.pca$rotation[, pc1], dat.pca$rotation[, pc2], xlab = paste0("PC", pc1), ylab = paste0("PC", pc2),  pch = ".", cex = 2, main = paste("Processing method:", dat.type, "Tissue:", jtiss))
text(dat.pca$rotation[, pc1], dat.pca$rotation[, pc2], labels = paste("ZT", colnames(M)), col = cols)
x <- dat.pca$x[, 1]
head(x[order(x, decreasing = TRUE)], n = 50)
