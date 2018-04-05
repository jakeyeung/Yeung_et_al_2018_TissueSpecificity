# 2016-07-14
# Jake Yeung

rm(list=ls())

setwd("~/projects/tissue-specificity")

library(ggplot2)
library(ggrepel)
library(dplyr)
library(hash)
library(reshape2)
library(hash)
library(PMA)

source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/PlotFunctions.R")
source("scripts/functions/SitecountsFunctions.R")
source("scripts/functions/LdaFunctions.R")
source("scripts/functions/RemoveP2Name.R")
source("scripts/functions/GetTopMotifs.R")
source("scripts/functions/HandleMotifNames.R")
source("scripts/functions/ListFunctions.R")
source("scripts/functions/PlotActivitiesFunctions.R")
source("scripts/functions/LoadActivitiesLong.R")
source("scripts/functions/LiverKidneyFunctions.R")
source("scripts/functions/SvdFunctions.R")
source("scripts/functions/NcondsAnalysisFunctions.R")
source("scripts/functions/FisherTestSitecounts.R")
source("scripts/functions/ColorFunctions.R")


# Function ----------------------------------------------------------------

BinPhases <- function(p){
  if (p > 24 | p < 0) stop("Not coded for phase > 24 or phase < 0")
  if (p > 20 | p <= 2){  # use OR to over
    bin <- "20-2"
  } else if (p > 2 & p <= 8){
    bin <- "2-8"
  } else if (p > 8 & p <= 14){
    bin <- "8-14"
  } else if (p > 14 & p <= 20){
    bin <- "14-20"
  }
}

# Main --------------------------------------------------------------------



pldarobj <- "/home/yeung/projects/tissue-specificity/Robjs/liver_kidney_atger_nestle/plda_robjs/2D.posterior.multigene.distfilt.morenonliv.bugfixed.liverWTKO.40000.cutoff.3.cutofflow0.method.g=1001.Robj"
# pldarobj <- "/home/yeung/projects/tissue-specificity/Robjs/liver_kidney_atger_nestle/plda_robjs/2D.posterior.multigene.distfilt.morenonliv.bugfixed.liverWTKO.50000.cutoff.3.cutofflow0.method.g=1001.model.Liver_SV129.Robj"
# pldarobj <- "/home/yeung/projects/tissue-specificity/Robjs/liver_kidney_atger_nestle/plda_robjs/2D.posterior.multigene.distfilt.morenonliv.bugfixed.liverWTKO.40000.cutoff.3.cutofflow0.5.method.g=1001.Robj"
load(pldarobj)


load("/home/yeung/projects/tissue-specificity/Robjs/liver_kidney_atger_nestle/fits.long.multimethod.filtbest.staggeredtimepts.bugfixed.Robj")
fits.long.filt.orig <- fits.long.filt
jmeth <- "g=1001"
fits.long.filt <- subset(fits.long.filt, method == jmeth)

# compare with Hogenesch
load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.15.phase_sd_maxdiff_avg.Robj")  # fits.best
load("Robjs/dat.long.fixed_rik_genes.Robj")
dat.long.hog <- dat.long
load("Robjs/liver_kidney_atger_nestle/dat.long.liverkidneyWTKO.bugfixed.Robj")
load("Robjs/fits.relamp.Robj")
load("Robjs/liver_kidney_atger_nestle/fits.bytiss.Robj")
fits.bytiss <- subset(fits.bytiss, !is.na(gene))
dat.orig <- dat.long
dat.long <- CollapseTissueGeno(dat.long, keep.tissue.col = TRUE)
dat.long <- StaggeredTimepointsLivKid(dat.long)

load("Robjs/liver_kidney_atger_nestle/dat.freq.bugfixed.Robj")

jsub.hog <- subset(fits.relamp, tissue == "Liver")
phase.hog <- hash(as.character(jsub.hog$gene), jsub.hog$phase)
jsub.atg <- subset(fits.bytiss, tissue == "Liver_SV129") 
phase.atg <- hash(as.character(jsub.atg$gene), jsub.atg$phase)
amp.atg <- hash(as.character(jsub.atg$gene), jsub.atg$amp)
pval.atg <- hash(as.character(jsub.atg$gene), jsub.atg$pval)



# Set up ------------------------------------------------------------------

mat.fgbg.lab.lst.3 <- SetUpMatForLda(mat.fg, mat.bgnonliver, mat.bg, has.peaks = TRUE)
mat.fgbg.3 <- mat.fgbg.lab.lst.3$mat.fgbg; labels3 <- mat.fgbg.lab.lst.3$labels
colnames(mat.fgbg.3) <- sapply(colnames(mat.fgbg.3), RemoveP2Name)
colnames(mat.fgbg.3) <- sapply(colnames(mat.fgbg.3), function(cname){
  return(RemoveCommasBraces(cname))
}, USE.NAMES = FALSE)


rhyth.motifs <- sapply(GetTopMotifs("rhythmic"), RemoveP2Name, USE.NAMES = FALSE)
rhyth.motifs <- c(rhyth.motifs, c("SRF"))
rhyth.motifs <- rhyth.motifs[which( ! rhyth.motifs %in% c("ATF6"))]
tissue.motifs <- sapply(GetTopMotifs("tissue"), RemoveP2Name, USE.NAMES = FALSE)
tissue.motifs <- c(tissue.motifs, c("ATF5_CREB3", "ATF6"))
tissue.motifs <- tissue.motifs[which( ! tissue.motifs %in% c("SRF"))]

# remove tissue motifs in rhyth
rhyth.motifs <- rhyth.motifs[which(!rhyth.motifs %in% intersect(rhyth.motifs, tissue.motifs))]

# cross prods
mat.rhyth3 <- subset(mat.fgbg.3, select = intersect(rhyth.motifs, colnames(mat.fgbg.3)))
mat.tiss3 <- subset(mat.fgbg.3, select = intersect(tissue.motifs, colnames(mat.fgbg.3)))
mat.rhythtiss3 <- CrossProductTwoSets(mat.rhyth3, mat.tiss3)

mat.fgbg.cross.rhythtiss3 <- cbind(mat.fgbg.3, mat.rhythtiss3)
# remove columns with 0 variance
mat.fgbg.cross.rhythtiss3[which(colSums(mat.fgbg.cross.rhythtiss3) == 0)] <- list(NULL)

jlambda <- 0.025  # liv only
jlambda <- 0.035  # liv only
out.cross.rhythtiss3 <- PenalizedLDA(mat.fgbg.cross.rhythtiss3, labels3, lambda = jlambda, K = 2, standardized = FALSE)


# Plot pretty -------------------------------------------------------------


# plot pretty
vec.length <- sqrt(out.cross.rhythtiss3$discrim[, 1]^2 + out.cross.rhythtiss3$discrim[, 2]^2)

jsize.cutoff <- 0.1
jsize.pairs.cut <- sapply(vec.length, function(jsize){
  if (jsize > jsize.cutoff){
    return(jsize)
  } else {
    return(0)
  }
})

labels <- names(out.cross.rhythtiss3$x)
labels.cut <- mapply(function(jlab, jsize){
  if (jsize <= 0){
    return("")
  } else {
    return(jlab)
  }
}, labels, jsize.pairs.cut)

dat.plot <- data.frame(x = out.cross.rhythtiss3$discrim[, 1],
                       y = out.cross.rhythtiss3$discrim[, 2],
                       motif = labels.cut,
                       vec.length = vec.length,
                       vec.length.cut = jsize.pairs.cut)
dat.labs <- subset(dat.plot, vec.length.cut > 0)

m <- ggplot(dat.plot, aes(x = x, y = y)) + 
  geom_point(size = 0.01) + 
  geom_text_repel(data = dat.labs, aes(x = x, y = y, label = motif), size = 2.5) + 
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + theme_bw() + 
  theme(aspect.ratio = 0.33, legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  xlab("Motif loadings separating liver and kidney DHS peaks") + ylab("Motif loadings separating rhythmic and flat DHS peaks")
print(m)


# Downstream --------------------------------------------------------------




jlabs <-  c("Liv rhyth", "Nonrhyth", "Liv flat")

jlabs.genes <-  c("#genes=", "#genes=", "#genes=")
counts.genes <- c(length(unique(mat.fg$gene)), length(unique(mat.bgnonliver$gene)), length(unique(mat.bg$gene)))
jlabs.counts <- paste(jlabs, counts.genes)

jlabs.peaks <-  c("Liv rhyth #peaks=", "Nonrhyth #peaks=", "Liv flat #peaks=")
counts.peaks <- c(length(mat.fg$gene), length(mat.bgnonliver$gene), length(mat.bg$gene))
jlabs.pcounts <- paste(jlabs.peaks, counts.peaks)

labs.annot <- paste(jlabs.counts, jlabs.pcounts, sep = "\n")

# count motifs
jlabs.counts <- paste(jlabs, "#peaks=", table(labels3))
mat.tmp <- mat.fgbg.cross.rhythtiss3
mat.tmp$peakgene <- rownames(mat.tmp)
jlabs.vec <- sapply(labels3, function(l){
  jlabs.counts[[l]]
})
mat.tmp$gene <- sapply(mat.tmp$peakgene, function(pg) strsplit(pg, ";")[[1]][[2]])
mat.tmp$label <- jlabs.vec
mat.tmp$label.n <- labels3
mat.tmp$phase <- sapply(as.character(mat.tmp$gene), function(g) ifelse(!is.null(phase.atg[[g]]), phase.atg[[g]], NA))
mat.tmp$amp <- sapply(as.character(mat.tmp$gene), function(g) ifelse(!is.null(amp.atg[[g]]), amp.atg[[g]], NA))
mat.tmp$pval <- sapply(as.character(mat.tmp$gene), function(g) ifelse(!is.null(pval.atg[[g]]), pval.atg[[g]], NA))

# add phase to label
mat.tmp$genephase <- paste(mat.tmp$gene, signif(mat.tmp$phase, 2), sep = ";")
# mat.tmp$genephaseamp <- paste()
mat.tmp$peakgenephase <- paste(mat.tmp$peak, signif(mat.tmp$phase, 2), sep = ";")

cnames <- c("label.n", "label", "peakgene", "gene", "phase", "amp", "pval", "genephase", "peakgenephase")
mat.long <- melt(mat.tmp, id.vars = cnames, variable.name = "motif", value.name = "sitecount")

# assign color to gene based on amp, phase, pval
mat.long$col <- PhaseAmpPvalToColor(mat.long$phase, mat.long$amp, mat.long$pval, rotate.hr = -8, pval.k = 0.5)

jmotif <- "FOXA2"
jmotif <- "RORA;ATF5_CREB3"
jmotif <- "RORA;FOXA2"
jmotif <- "CEBPA.B_DDIT3"
jmotif <- "bHLH_family;ATF5_CREB3"
jmotif <- "RORA;CUX2"
txtsize <- 4

# label outliers http://stackoverflow.com/questions/33524669/labeling-outliers-of-boxplots-in-r
subset(mat.long, motif == jmotif & sitecount > 0) %>%
  group_by(label) %>%
  mutate(outlier = ifelse(is_top_N(sitecount, N = 10), peakgenephase, as.numeric(NA))) %>%
  ggplot(., aes(x = label, y = sitecount, colour = col)) +
  geom_boxplot(aes(colour = NULL), outlier.shape = NA) +
  # geom_text(aes(label = outlier), na.rm = TRUE, hjust = -0.3) + 
  geom_text_repel(aes(label = outlier), size = txtsize, na.rm = TRUE) +
  geom_point() + 
  # geom_text_repel(aes(label = outlier), na.rm = TRUE) + 
  ggtitle(jmotif) + 
  theme_bw() + 
  xlab("") + ylab("Sitecount") + 
  scale_colour_identity()


# Remove redundancies -----------------------------------------------------

fits.liv <- subset(fits.bytiss, tissue == "Liver_SV129")
amp.hash <- hash(as.character(fits.liv$gene), fits.liv$amp)
phase.hash <- hash(as.character(fits.liv$gene), fits.liv$phase)
pval.hash <- hash(as.character(fits.liv$gene), fits.liv$pval)

discrim <- out.cross.rhythtiss3$discrim
rownames(discrim) <- names(out.cross.rhythtiss3$x)
discrim.hits <- discrim[which(apply(discrim, 1, min) > 0.0125), ]
motifs.hits <- rownames(discrim.hits)


mat.M <- mat.fgbg.cross.rhythtiss3[which(labels3 == 1), motifs.hits]
# Plot loadings
mat.M <- mat.M[, grepl(";", colnames(mat.M))]
mat.M <- mat.M[, !grepl("FOX\\.F1\\.F2", colnames(mat.M))]
mat.M <- t(scale(t(mat.M), center = FALSE, scale = FALSE))
mat.M <- mat.M[which(rowSums(mat.M) > 0), ]

mat.pca <- prcomp(mat.M, center = FALSE, scale. = FALSE)
# biplot
biplot(mat.pca)
# separately
pc1 <- 1
pc2 <- 2
plot(loadings[, pc1], loadings[, pc2])
text(loadings[, pc1], loadings[, pc2], rownames(loadings))
plot(motif.space[, pc1], motif.space[, pc2])
text(motif.space[, pc1], motif.space[, pc2], rownames(motif.space))

# Find association by bins!
loadings <- mat.pca$rotation
rownames(loadings) <- motifs.hits
motif.space <- mat.pca$x
rownames(motif.space) <- rownames(mat.M)
# second component ranks by rhythmic motif content
motif.pc1 <- motif.space[, pc1]
motif.pc2 <- motif.space[, pc2]
# motif.pc2 <- sort(motif.pc2)
motif.pc2.df <- data.frame(pc1 = motif.pc1,
                           pc2 = motif.pc2, 
                           i = seq(length(motif.pc2)), 
                           genepeak = names(motif.pc2), 
                           gene = sapply(names(motif.pc2), function(m) strsplit(m, ";")[[1]][[2]]))
motif.pc2.df$amp <- sapply(as.character(motif.pc2.df$gene), function(g) amp.hash[[g]])
motif.pc2.df$phase <- sapply(as.character(motif.pc2.df$gene), function(g) phase.hash[[g]])
motif.pc2.df$pval <- sapply(as.character(motif.pc2.df$gene), function(g) pval.hash[[g]])

ggplot(motif.pc2.df, aes(x = amp, y = pc2)) + geom_point()
ggplot(motif.pc2.df, aes(x = phase, y = pc2)) + geom_point(alpha = 0.25)
ggplot(motif.pc2.df, aes(x = phase, y = pc1)) + geom_point(alpha = 0.25)

# bin into phases then do boxplot
motif.pc2.df$phase.bin <- sapply(motif.pc2.df$phase, BinPhases)

ggplot(subset(motif.pc2.df, abs(pc2) > 0.1), aes(x = phase.bin, y = pc2)) + geom_boxplot() 
ggplot(subset(motif.pc2.df, abs(pc2) > 0.02), aes(x = phase.bin, y = pc2)) + geom_boxplot() 

# ggplot(motif.pc2.df, aes(x = phase.bin, y = pc2)) + geom_violin() 
# ggplot(motif.pc2.df, aes(x = phase.bin, y = pc2)) + geom_boxplot() 
# ggplot(motif.pc2.df, aes(x = phase.bin, y = pc1)) + geom_violin() 
# ggplot(motif.pc2.df, aes(x = phase.bin, y = pc1)) + geom_boxplot() 

plot(seq(motif.pc2), motif.pc2)
text(seq(motif.pc2), motif.pc2, names(motif.pc2))

# Do for penalized
# mat.M.outlierremoved <- mat.M[!grepl("Pde9a", rownames(mat.M)), ]
mat.pmd <- PMD(mat.M, type = c("standard"), sumabs = 0.4, center=FALSE, rnames = rownames(mat.M), cnames = colnames(mat.M), K=2)
print(mat.pmd)

rownames(mat.pmd$u) <- mat.pmd$rnames
rownames(mat.pmd$v) <- mat.pmd$cnames

plot(x = mat.pmd$v[, 1], y = mat.pmd$v[, 2])
text(x = mat.pmd$v[, 1], y = mat.pmd$v[, 2], labels = mat.pmd$cnames)

plot(x = mat.pmd$u[, 1], y = mat.pmd$u[, 2])
text(x = mat.pmd$u[, 1], y = mat.pmd$u[, 2], labels = mat.pmd$rnames)

biplot(x = mat.pmd$u, y = mat.pmd$v, xlabs = mat.pmd$rnames, ylabs = mat.pmd$cnames)

motif.pmd.df <- data.frame(pc1 = mat.pmd$u[, 1],
                           pc2 = mat.pmd$u[, 2], 
                           genepeak = rownames(mat.pmd$u), 
                           gene = sapply(rownames(mat.pmd$u), function(m) strsplit(m, ";")[[1]][[2]]))
motif.pmd.df$amp <- sapply(as.character(motif.pmd.df$gene), function(g) amp.hash[[g]])
motif.pmd.df$phase <- sapply(as.character(motif.pmd.df$gene), function(g) phase.hash[[g]])
motif.pmd.df$pval <- sapply(as.character(motif.pmd.df$gene), function(g) pval.hash[[g]])

motif.pmd.df$phase.bin <- sapply(motif.pmd.df$phase, BinPhases)

ggplot(subset(motif.pmd.df, abs(pc2) > 0.1), aes(x = phase.bin, y = pc2)) + geom_boxplot() 
ggplot(subset(motif.pmd.df, abs(pc1) > 0), aes(x = phase.bin, y = pc1)) + geom_boxplot() 



# Why is Ept1 do not have proper tissue-specifi peaks? --------------------

if (!exists("S.long")) load("Robjs/S.long.multigene.filt.50000.Robj", v=T)
S.long <- subset(S.long, tissue %in% c("Kidney", "Liver"))

# Heatmap -----------------------------------------------------------------


# correlate 
library(gplots)
source("scripts/functions/PlotFunctions.R")

jtitle <- ""
blackend <- 0.25
minval <- 0
maxval <- 1
my.palette <- colorRampPalette(c("black", "yellow"))(n = 300)
# # (optional) defines the color breaks manually for a "skewed" color transition
col.breaks = c(seq(minval, blackend, length=150),  # black
               seq(blackend + 0.01, maxval, length=151))  # yellow
jdendro <- "both"

heatmap.2(as.matrix(mat.M),
          col=my.palette,
          breaks = col.breaks,
          scale="none",
          key=T,
          keysize=1.5,
          density.info = "density",
          trace="none",
          cexCol=1.4,
          labRow=NA,
          Rowv=TRUE,
          Colv=TRUE,
          dendrogram = jdendro,
          main = jtitle,
          margins=c(15, 5))
