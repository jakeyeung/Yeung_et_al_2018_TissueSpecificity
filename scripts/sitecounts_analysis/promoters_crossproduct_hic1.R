# 2016-04-012
# promoters_crossproduct_hic1.R

library(hash)
library(ggplot2)
library(dplyr)

# Sources -----------------------------------------------------------------

source('/home/yeung/projects/ridge-regression/ridgeInR.R')
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/GetClockGenes.R")
source("scripts/functions/LoadActivitiesPlotSvd.R")
source("scripts/functions/LdaFunctions.R")
load("Robjs/dat.long.fixed_rik_genes.Robj", v=T)
load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.15.phase_sd_maxdiff_avg.Robj", v=T)

cgenes <- GetClockGenes()
# Infs --------------------------------------------------------------------

site <- "/home/yeung/projects/tissue-specificity/data/sitecounts/motevo/sitecount_matrix_geneids"
N <- read.table.handlerows(site)
N.geneid <- read.table(site, header=TRUE)
genes <- N.geneid$Gene.ID

# Observe HIC1 ------------------------------------------------------------

jmotif <- "HIC1.p2"
plot(density(N[[jmotif]]))

N.sub <- subset(N, select = c("Gene.ID", jmotif))

N.sub <- N.sub[order(N.sub$HIC1.p2, decreasing = TRUE), ]

hic1.hash <- hash(as.character(N.sub$Gene.ID), as.character(N.sub$HIC1.p2))  # vectorize probably faster?

# Do cross prods ----------------------------------------------------------

hic1.vec <- N$HIC1.p2
N.cross <- sweep(x = N, MARGIN = 1, STATS = hic1.vec, FUN = "*")
N.cross$HIC1.p2 <- NULL  # ignore self

# rename
colnames(N.cross) <- sapply(colnames(N.cross), function(cname) paste0(cname, ";HIC1.p2"))

# check top crosses in core clocks?
N.cross.cgenes <- N.cross[which(genes %in% cgenes), ]
apply(N.cross.cgenes, 1, max)

# RBIND
N.withcross <- cbind(N.geneid, N.cross)
N.crossonly <- cbind(N.geneid$Gene.ID, N.cross)

write.table(N.withcross, file = "/home/yeung/projects/tissue-specificity/data/sitecounts/motevo/sitecount_matrix_geneids.HIC1_cross.mat", quote = FALSE, row.names = FALSE, sep = "\t")
write.table(N.crossonly, file = "/home/yeung/projects/tissue-specificity/data/sitecounts/motevo/sitecount_matrix_geneids.HIC1_cross_only.mat", quote = FALSE, row.names = FALSE, sep = "\t")


# Run MARA from script ----------------------------------------------------



# Analyze MARA output -----------------------------------------------------

LoadActivitiesPlotSvd("/home/yeung/projects/tissue-specificity/results/MARA/bic_modules/TissueWide.centeredTRUE", 1)
LoadActivitiesPlotSvd("/home/yeung/projects/tissue-specificity/results/MARA/bic_modules_crossprod/TissueWide.centeredTRUE", 1)
LoadActivitiesPlotSvd("/home/yeung/projects/tissue-specificity/results/MARA/all_genes_crossprod_singlefacs_and_hic1/expressed_genes_deseq_int.centeredTRUE", 2)

LoadActivitiesPlotSvd(inf <- "/home/yeung/projects/tissue-specificity/results/MARA/all_genes_crossprod_hic1/expressed_genes_deseq_int.centeredTRUE", 1)


# Does tissue wide genes correlate with HIC1 ------------------------------

jsub <- subset(fits.best, n.rhyth >= 8)
# jsub <- subset(fits.best, model == "Adr")
genes.tw <- as.character(unique(jsub$gene))
# assign avg amp
genes.tw.amp <- hash(as.character(jsub$gene), jsub$amp.avg)

N.sub <- subset(N.geneid, Gene.ID %in% genes.tw)
df.hic1.amp <- data.frame(gene=N.sub$Gene.ID, hic1.sitecount=N.sub$HIC1, avg.amp=sapply(as.character(N.sub$Gene.ID), function(g) genes.tw.amp[[g]]))

ggplot(df.hic1.amp, aes(x=hic1.sitecount, y=avg.amp, label=gene)) + geom_text()

# cross product with clock factors
clock.motifs <- c("RORA.p2", "HSF1.2.p2", "bHLH_family.p2", "NFIL3.p2", "SRF.p3")
clock.motifs <- c("RORA.p2", "bHLH_family.p2", "NFIL3.p2")
# clock.motifs <- c("HNF4A_NR2F1.2.p2")

N.sub.clocks <- subset(N.sub, select = clock.motifs)
HIC1.vec <- N.sub$HIC1.p2
gene.vec <- N.sub$Gene.ID

# cross only clock motifs only
N.sub.clocks.cross <- sweep(x = N.sub.clocks, MARGIN = 1, STATS = HIC1.vec, FUN = "*")
# N.sub.clocks.cross <- sweep(x = N.sub.clocks, MARGIN = 1, STATS = 1, FUN = "*")

N.sub.clocks.cross.max <- apply(N.sub.clocks.cross, 1, max)
df.hic1.clockcross.amp <- data.frame(gene=N.sub$Gene.ID, hic1.sitecount=N.sub.clocks.cross.max, avg.amp=sapply(as.character(N.sub$Gene.ID), function(g) genes.tw.amp[[g]]))
ggplot(df.hic1.clockcross.amp, aes(x=hic1.sitecount, y=avg.amp, label=gene)) + geom_text() + expand_limits(y=0)
ggplot(df.hic1.clockcross.amp, aes(x=hic1.sitecount, y=avg.amp, label=gene)) + geom_point(alpha=0.2) + expand_limits(y=0)


# Do cross prdoucts for all pairs -----------------------------------------

N.sub.cross.all <- CrossProduct(mat = subset(N.sub, select = -Gene.ID))
N.sub.cross.all <- cbind(N.sub, N.sub.cross.all)  # add single facs

write.table(N.sub.cross.all, file = "/home/yeung/projects/tissue-specificity/data/sitecounts/motevo/sitecount_matrix_geneids.cross_all.tissuewide.mat", quote = FALSE, row.names = FALSE, sep = "\t")


# Load MARA results -------------------------------------------------------

LoadActivitiesPlotSvd("/home/yeung/projects/tissue-specificity/results/MARA/bic_modules_crossprod_all/TissueWide.centeredTRUE", comp = 1)


# Do cross products on a list of motifs -----------------------------------

# cross.motifs <- c("RORA.p2", "HSF1.2.p2", "bHLH_family.p2", "NFIL3.p2", "SRF.p3", "HIC1.p2", "TFDP1.p2", "ATF2.p2", "MYBL2.p2", "SP1.p2", "PAX5.p2", "NFIX.p2", "EP300.p2", "NR3C1.p2", "HLF.p2", "NFE2L2.p2")
cross.motifs <- c("RORA.p2", "NFIL3.p2", "bHLH_family.p2", "HSF1.2.p2")

N.crosses <- lapply(cross.motifs, function(jmotif, N){
  cross.vec <- N[[jmotif]]
  N.cross <- sweep(x = N, MARGIN = 1, STATS = cross.vec, FUN = "*")
  N.cross[[jmotif]] <- NULL
  colnames(N.cross) <- sapply(colnames(N.cross), function(cname) paste0(cname, ";", jmotif))
  return(N.cross)
}, N)
N.crosses <- do.call(cbind, N.crosses)

# add single factors
N.crosses <- cbind(N.geneid, N.crosses)

write.table(N.crosses, file = "/home/yeung/projects/tissue-specificity/data/sitecounts/motevo/sitecount_matrix_geneids.cross_clockfactors.mat", quote = FALSE, row.names = FALSE, sep = "\t")


# Load MARA ---------------------------------------------------------------

# LoadActivitiesPlotSvd("/home/yeung/projects/tissue-specificity/results/MARA/bic_modules_crossprod_several/TissueWide.centeredTRUE", comp = 1)
LoadActivitiesPlotSvd("/home/yeung/projects/tissue-specificity/results/MARA/bic_modules_crossprod_rora/TissueWide.centeredTRUE", comp = 1)
LoadActivitiesPlotSvd("/home/yeung/projects/tissue-specificity/results/MARA/bic_modules_crossprod_clockfactors/TissueWide.centeredTRUE", comp=1)

