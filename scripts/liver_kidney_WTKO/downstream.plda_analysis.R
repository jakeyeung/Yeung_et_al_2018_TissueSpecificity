# 2016-07-14
# Jake Yeung

rm(list=ls())

library(ggplot2)
library(ggrepel)
library(dplyr)
library(hash)
library(reshape2)

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
source("scripts/functions/PlotUCSC.R")
source("scripts/functions/LiverKidneyFunctions.R")


# Functions ---------------------------------------------------------------



# Load --------------------------------------------------------------------

# LIVER
# load("/home/yeung/projects/tissue-specificity/Robjs/liver_kidney_atger_nestle/plda_robjs/2D.posterior.multigene.distfilt.morenonliv.bugfixed.liverWTKO.40000.cutoff.3.cutofflow0.method.g=1001.Robj", v=T)
# load("/home/yeung/projects/tissue-specificity/Robjs/liver_kidney_atger_nestle/plda_robjs/2D.posterior.multigene.distfilt.morenonliv.bugfixed.liverWTKO.40000.cutoff.3.cutofflow0.method.BIC.Robj", v=T)
pldarobj <- "/home/yeung/projects/tissue-specificity/Robjs/liver_kidney_atger_nestle/plda_robjs/2D.posterior.multigene.distfilt.morenonliv.bugfixed.liverWTKO.40000.cutoff.3.cutofflow0.method.g=1001.Robj"
jmodel <- "Liver_SV129"

# KIDNEY
# pldarobj <- "/home/yeung/projects/tissue-specificity/Robjs/liver_kidney_atger_nestle/plda_robjs/2D.posterior.multigene.distfilt.morenonliv.bugfixed.liverWTKO.40000.cutoff.3.cutofflow0.method.g=1001.model.Kidney_SV129.Robj"
# pldarobj <- "/home/yeung/projects/tissue-specificity/Robjs/liver_kidney_atger_nestle/plda_robjs/2D.posterior.multigene.distfilt.morenonliv.bugfixed.liverWTKO.40000.cutoff.2.5.cutofflow0.5.method.g=1001.model.Kidney_SV129.Robj"
# jmodel <- "Kidney_SV129"


load(pldarobj, v=T)

load("/home/yeung/projects/tissue-specificity/Robjs/liver_kidney_atger_nestle/fits.long.multimethod.filtbest.staggeredtimepts.bugfixed.Robj", v=T)

# compare with Hogenesch
load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.15.phase_sd_maxdiff_avg.Robj", v=T)  # fits.best
load("Robjs/dat.long.fixed_rik_genes.Robj", v=T)
dat.long.hog <- dat.long
load("Robjs/liver_kidney_atger_nestle/dat.long.liverkidneyWTKO.bugfixed.Robj", v=T)
load("Robjs/fits.relamp.Robj", v=T)
load("Robjs/liver_kidney_atger_nestle/fits.bytiss.Robj", v=T)
dat.orig <- dat.long
dat.long <- CollapseTissueGeno(dat.long, keep.tissue.col = TRUE)
dat.long <- StaggeredTimepointsLivKid(dat.long)


liver.genes <- as.character(subset(fits.long.filt, method == "g=1001" & model == "Liver_SV129")$gene)
liver.genes.hog <- as.character(subset(fits.best, model == "Liver")$gene)

fits.livhog <- subset(fits.best, gene %in% liver.genes)
jsub.hog <- subset(fits.relamp, tissue == "Liver" & gene %in% liver.genes)
phase.hog <- hash(as.character(jsub.hog$gene), jsub.hog$phase)
jsub.atg <- subset(fits.bytiss, tissue == "Liver_SV129" & gene %in% liver.genes) 
phase.atg <- hash(as.character(jsub.atg$gene), jsub.atg$phase)

phase.shift <- data.frame(gene = liver.genes)
phase.shift$phase.hog <- sapply(as.character(phase.shift$gene), function(g){
  p <- phase.hog[[g]]
  if (is.null(p)){
    p <- NA
  }
  return(p)
})
phase.shift$phase.atg <- sapply(as.character(phase.shift$gene), function(g){
  p <- phase.atg[[g]]
  if (is.null(p)){
    p <- NA
  }
  return(p)
})

ggplot(phase.shift[complete.cases(phase.shift), ], aes(x = phase.hog, y = phase.atg, label = gene)) + 
  geom_point(alpha = 0.25) + 
  xlim(c(0, 24)) + ylim(c(0, 24)) + theme_bw() + xlab("Phase in DD, AL") + ylab("Phase in LD, RF")


# Find motif enrichment ---------------------------------------------------


# colnames(mat.fg) <- sapply(colnames(mat.fg), RemoveP2Name)
# colnames(mat.fg) <- sapply(colnames(mat.fg), RemoveCommasBraces)
# colnames(mat.bgnonliver) <- sapply(colnames(mat.bgnonliver), RemoveP2Name)
# colnames(mat.bgnonliver) <- sapply(colnames(mat.bgnonliver), RemoveCommasBraces)
# colnames(mat.bg) <- sapply(colnames(mat.bg), RemoveP2Name)
# colnames(mat.bg) <- sapply(colnames(mat.bg), RemoveCommasBraces)

# remove a gene to see what happens
# mat.fg <- subset(mat.fg, gene != "Insig2")

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
out.cross.rhythtiss3 <- PenalizedLDA(mat.fgbg.cross.rhythtiss3, labels3, lambda = jlambda, K = 2, standardized = FALSE)

# plot pretty
vec.length <- sqrt(out.cross.rhythtiss3$discrim[, 1]^2 + out.cross.rhythtiss3$discrim[, 2]^2)
jsize.pairs <- vec.length * 5 + 0.01
plot(out.cross.rhythtiss3$discrim[, 1], out.cross.rhythtiss3$discrim[, 2], pch = ".",
     xlab = "Liver-specific DHS vs Kidney-specific DHS", ylab = "Liver-specific rhythmic DHS vs Nonrhythmic DHS",
     main="Motifs that separate between tissues (x-axis) and rhythmic in liver (y-axis)")
text(out.cross.rhythtiss3$discrim[, 1], out.cross.rhythtiss3$discrim[, 2], names(out.cross.rhythtiss3$x), cex = jsize.pairs)
abline(v = 0); abline(h = 0)

jsize.cutoff <- 0.09
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
  theme(aspect.ratio = 0.25, legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  xlab("Motif loadings separating liver and kidney DHS peaks") + ylab("Motif loadings separating rhythmic and flat DHS peaks")
print(m)

jlabs <-  c("Liv rhyth", "Nonrhyth", "Liv flat")
gcounts <- c(length(unique(mat.fg$gene)), length(unique(mat.bgnonliver$gene)), length(unique(mat.bg$gene)))

BoxplotLdaOut(out.cross.rhythtiss3, jdim = 2, horizontal = FALSE, axis.names = jlabs, jtitle = "Yaxis: Liv rhyth vs others")
BoxplotLdaOut(out.cross.rhythtiss3, jdim = 1, horizontal = TRUE, axis.names = jlabs, jtitle = "Xaxis: liver vs nonliver")

print(paste("Number of genes in mat.fg:", length(unique(mat.fg$gene))))
print(paste("Number of genes in mat.bgnonliver:", length(unique(mat.bgnonliver$gene))))
print(paste("Number of genes in mat.bg:", length(unique(mat.bg$gene))))


# Find genes with motifs enriched  ----------------------------------------------------

# count motifs
jlabs.counts <- paste(jlabs, "n=", table(labels3))
mat.tmp <- mat.fgbg.cross.rhythtiss3
mat.tmp$peakgene <- rownames(mat.tmp)
jlabs.vec <- sapply(labels3, function(l){
  jlabs.counts[[l]]
})
mat.tmp$gene <- sapply(mat.tmp$peakgene, function(pg) strsplit(pg, ";")[[1]][[2]])
mat.tmp$label <- jlabs.vec
mat.tmp$label.n <- labels3
mat.tmp$phase <- sapply(as.character(mat.tmp$gene), function(g) ifelse(!is.null(phase.atg[[g]]), phase.atg[[g]], NA))

mat.long <- melt(mat.tmp, id.vars = c("label.n", "label", "peakgene", "gene", "phase"), variable.name = "motif", value.name = "sitecount")

jmotif <- c("NR4A2;CUX2")
jmotif <- c("ESRRA")
jmotif <- c("FOX.F1.F2.J1.;CUX2")
jmotif <- c("CEBPA.B_DDIT3")
jmotif <- c("SRF;CUX2")
jmotif <- c("FOXA2")
jmotif <- c("PATZ1")
jmotif <- c("RORA;HNF4A_NR2F1.2")
jmotif <- c("RORA;GATA6")
jmotif <- c("GATA6")
jmotif <- c("XBP1;ONECUT1.2")
jmotif <- c("RORA;FOXA2")
jmotif <- c("RORA;ONECUT1.2")
jmotif <- c("RORA;CUX2")
jmotif <- c("NR4A2;ONECUT1.2")
jmotif <- c("TFCP2")
jmotif <- c("RORA;ATF5_CREB3")
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# label outliers http://stackoverflow.com/questions/33524669/labeling-outliers-of-boxplots-in-r
subset(mat.long, motif == jmotif) %>%
  group_by(label) %>%
  mutate(outlier = ifelse(is_top_N(sitecount, N = 10), peakgene, as.numeric(NA))) %>%
  ggplot(., aes(x = label, y = sitecount)) +
  geom_boxplot() +
  # geom_text(aes(label = outlier), na.rm = TRUE, hjust = -0.3) + 
  geom_text_repel(aes(label = outlier), na.rm = TRUE) + 
  ggtitle(jmotif) + 
  theme_bw() + 
  xlab("") + ylab("Sitecount")

# show phase relationship
ggplot(subset(mat.long, motif == jmotif), aes(y = sitecount, x = phase)) + geom_point() + ggtitle(jmotif) + theme_bw()

# faster
ggplot(subset(mat.long, motif == jmotif), aes(y = sitecount, x = label)) + geom_boxplot()

# densities hard to interpret
ggplot(subset(mat.long, motif == jmotif & label != "Nonrhyth"), aes(x = sitecount, fill = label)) + geom_density(alpha = 0.33) + scale_fill_manual(values = cbPalette)
ggplot(subset(mat.long, motif == jmotif), aes(x = sitecount, fill = label)) + geom_density(alpha = 0.33) + scale_fill_manual(values = cbPalette)

fg.sub <- mat.fgbg.cross.rhythtiss3[which(labels3 == 1), jmotif]
bg.sub <- mat.fgbg.cross.rhythtiss3[which(labels3 == 3), jmotif]

