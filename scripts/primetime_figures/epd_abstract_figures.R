# 2016-08-07
# Jake Yeung
# Do figures for EPD conference paper

rm(list=ls())

library(ggplot2)
library(PMA)
# detach("plyr", unload=TRUE)
library(dplyr)

setwd("/home/yeung/projects/tissue-specificity")

source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/NcondsFunctions.R")
source("scripts/functions/SvdFunctions.R")
source("scripts/functions/LoadActivitiesLong.R")
source("scripts/functions/LiverKidneyFunctions.R")
source("scripts/functions/PlotActivitiesFunctions.R")
source("scripts/functions/AnalyzeGeneEnrichment.R")
source("scripts/functions/FourierFunctions.R")
source("scripts/functions/GetTFs.R")
source("scripts/functions/LdaFunctions.R")
source("scripts/functions/HandleMotifNames.R")
source("scripts/functions/RemoveP2Name.R")
source("scripts/functions/GetTopMotifs.R")

remove.wfat <- TRUE
plot.dir <- "/home/yeung/projects/tissue-specificity/plots/primetime_plots_epd_abstract"

# Tissue-wide modules -----------------------------------------------------

tfs <- GetTFs(get.mat.only = TRUE)

jmeth <- "g=1001"
load("Robjs/liver_kidney_atger_nestle/dat.freq.bugfixed.Robj", v=T)  # LivKid
load("Robjs/liver_kidney_atger_nestle/fits.long.multimethod.filtbest.staggeredtimepts.bugfixed.annotated.Robj", v=T)
fits.long.filt <- subset(fits.long.filt, method == jmeth)
load("Robjs/liver_kidney_atger_nestle/dat.long.liverkidneyWTKO.bugfixed.Robj", v=T); dat.wtko <- dat.long; rm(dat.long)
load("Robjs/liver_kidney_atger_nestle/fits.bytiss.bugfixed.Robj", v=T)

load("Robjs/dat.long.fixed_rik_genes.Robj", v=T)  # hogenesch
# model selection with g1000
load("Robjs/nconds_g1000_11_tissues/fits_long.11_tiss_3_max.g1000.bestmodel.filteramp.0.15.Robj", v=T)
load("Robjs/dat.complex.fixed_rik_genes.Robj", v=T)
if (remove.wfat){
  dat.complex <- subset(dat.complex, tissue != "WFAT")
  dat.long <- subset(dat.long, tissue != "WFAT")
}



dat.wtko <- StaggeredTimepointsLivKid(dat.wtko)
dat.wtko.collapsed <- CollapseTissueGeno(dat.wtko)

# get input genes
genes.tw <- as.character(subset(fits.long, n.rhyth >= 8)$gene)
genes.tw.wtko <- as.character(subset(fits.long.filt, model %in% c("Liver_SV129,Kidney_SV129"))$gene)

# get regulators: hogenesch 
outbase <- "/home/yeung/projects/tissue-specificity/results/MARA.hogenesch"
outmain <- file.path(outbase, paste0("promoters.tissuewide.filteramp.0.15.mat"))
indir <- file.path(outmain, "expressed_genes_deseq_int.centeredTRUE")
act.long <- LoadActivitiesLong(indir, shorten.motif.name = TRUE)
# rename motifs based on the motifs with non-zero entries
omega <- 2 * pi / 24
act.complex <- subset(act.long, tissue != "WFAT") %>%
  group_by(gene, tissue) %>%
  do(ProjectToFrequency2(., omega, add.tissue=TRUE))
s.act <- SvdOnComplex(act.complex, value.var = "exprs.transformed")

# get regulators: liver kidney WT
indir <- "/home/yeung/projects/tissue-specificity/results/MARA.liver_kidney/promoters.Liver_SV129,Kidney_SV129.g=1001/atger_with_kidney.bugfixed"
act.wtko <- LoadActivitiesLongKidneyLiver(indir, shorten.motif.name = TRUE)
# rename motifs based on the motifs with non-zero entries
omega <- 2 * pi / 24
act.complex.wtko <- act.wtko %>%
  group_by(gene, tissue) %>%
  do(ProjectToFrequency2(., omega, add.tissue=TRUE))
s.act.wtko <- SvdOnComplex(act.complex.wtko, value.var = "exprs.transformed")

# plot 
comp <- 1
pdf(file.path(plot.dir, "tissue_wide_genes_and_regulators.pdf"))
  s.tw <- SvdOnComplex(subset(dat.complex, gene %in% genes.tw), value.var = "exprs.transformed")
  eigens.tw <- GetEigens(s.tw, period = 24, comp = comp, label.n = 15, eigenval = TRUE, adj.mag = TRUE, constant.amp = 4, peak.to.trough = TRUE)
  jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)
  print(eigens.tw$u.plot)
  multiplot(eigens.tw$u.plot, eigens.tw$v.plot, layout = jlayout)
  
  eigens.act <- GetEigens(s.act, period = 24, comp = comp, adj.mag = TRUE, constant.amp = 4, label.n = 25, jtitle = "", peak.to.trough = TRUE)
  print(eigens.act$u.plot)
  multiplot(eigens.act$u.plot, eigens.act$v.plot, cols = 2)
  
  s <- SvdOnComplex(subset(dat.freq, gene %in% genes.tw.wtko), value.var = "exprs.transformed")
  eigens <- GetEigens(s, period = 24, comp = comp, label.n = 15, eigenval = TRUE, adj.mag = TRUE, constant.amp = 4, peak.to.trough = TRUE)
  jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)
  multiplot(eigens$u.plot, eigens$v.plot, layout = jlayout)
  
  eigens.act <- GetEigens(s.act.wtko, period = 24, comp = comp, adj.mag = TRUE, constant.amp = 4, label.n = 25, jtitle = "", peak.to.trough = TRUE)
  print(eigens.act$u.plot)
  multiplot(eigens.act$u.plot, eigens.act$v.plot, cols = 2)
dev.off()


# Circadian flucutations in cellular heterogeneity ------------------------


exprs.thres <- 5
show.top.n <- 20
# ggplot(dat.mean.rnaseq, aes(x = exprs.mean)) + geom_density() + facet_wrap(~tissue)  # to find threshold

filt.tiss <- "WFAT"

indir <- "/home/yeung/projects/tissue-specificity/results/MARA/MARA_motevo_with_se.redo/activities"
act.long.hog.genomewide <- LoadActivitiesLong(indir, shorten.motif.name = TRUE)

dat.mean.rnaseq <- subset(dat.long, experiment == "rnaseq") %>%
  group_by(gene, tissue, experiment) %>%
  summarise(exprs.mean = mean(exprs))

pdf(file.path(plot.dir, "tissue_specific_hogenesch.pdf"))
  jtiss <- "BFAT"
  BFAT.genes <- as.character(subset(fits.long, model == jtiss)$gene)
  m <- PlotOverlayTimeSeries(dat.long, BFAT.genes, tissues = jtiss, jalpha = 0.05, jtitle = paste0(jtiss, "-specific rhythmic genes"))
  print(m)
  m2 <- PlotOverlayTimeSeries(subset(dat.long, experiment == "array"), BFAT.genes, tissues = jtiss, jalpha = 0.05, jtitle = paste0(jtiss, "-specific rhythmic genes"))
  print(m2)
  PlotMeanExprsOfModel(dat.mean.rnaseq, BFAT.genes, jtiss, avg.method = "mean")
  
  jtiss <- "Mus"
  Mus.genes <- as.character(subset(fits.long, model == jtiss)$gene)
  m <- PlotOverlayTimeSeries(dat.long, Mus.genes, tissues = jtiss, jscale = T, jalpha = 0.05, jtitle = paste0(jtiss, "-specific rhythmic genes"))
  print(m)
  m2 <- PlotOverlayTimeSeries(subset(dat.long, experiment == "array"), Mus.genes, tissues = jtiss, jalpha = 0.05, jtitle = paste0(jtiss, "-specific rhythmic genes"))
  print(m2)
  PlotMeanExprsOfModel(dat.mean.rnaseq, Mus.genes, jtiss, avg.method = "median")
  
  act.long <- subset(act.long, !tissue %in% filt.tiss)
  jmotifs.lst <- list("BFAT"="MEF2.A.B.C.D.", "Mus"="SPIB")
  for (jtiss in names(jmotifs.lst)){
    jmotif <- jmotifs.lst[[jtiss]]
    for (jexp in c("array", "rnaseq")){
      m <- PlotActivitiesWithSE(subset(act.long.hog.genomewide, gene == jmotif & tissue == jtiss & experiment == jexp)) + theme_bw() + theme(aspect.ratio = 1, legend.position = "none")
      m1 <- PlotActivitiesWithSE(subset(act.long.hog.genomewide, gene == jmotif & experiment == jexp)) + theme_bw() + theme(aspect.ratio = 1, legend.position = "none")
      print(m)
      print(m1)
    }
  }
  
  # Do GO enrichment
  for (jtiss in c("BFAT", "Mus")){
    genes.bg <- as.character(subset(dat.mean.rnaseq, exprs.mean > exprs.thres & tissue == jtiss)$gene)
    genes.fg <- as.character(subset(fits.long, model == jtiss)$gene)
    enrichment <- AnalyzeGeneEnrichment(genes.bg, genes.fg)
    enrichment$minuslogpval <- -log10(as.numeric(enrichment$classicFisher))
    enrichment <- OrderDecreasing(enrichment, jfactor = "Term", jval = "minuslogpval")
    enrichment <- enrichment[1:show.top.n, ]
    m1 <- ggplot(enrichment, aes(x = Term, y = minuslogpval)) + geom_bar(stat = "identity") + 
      ylab("-log10(P-value), Fisher's exact test") + 
      xlab("") +
      theme_bw() + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1), aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    print(m1)
  }
dev.off()

# Summarize modules -------------------------------------------------------

# filter top models
fits.count <- subset(fits.long.filt, method == jmeth & model != "") %>% group_by(model) %>% summarise(model.count = length(model))
fits.count <- fits.count[order(fits.count$model.count, decreasing = TRUE), ]
fits.count <- fits.count[1:10, ]
top.models <- as.character(fits.count$model)

fits.sub <- subset(fits.long.filt, model %in% top.models)
# reorder factor by fits.count
fits.sub$model <- factor(as.character(fits.sub$model), levels = top.models)


amp.thres <- seq(from = 0, to = max(fits.bytiss$amp), by = 0.15)
fits.counts.by.amp <- fits.sub %>%
  group_by(model) %>%
  do(NGenesByAmp.long(., amp.thres, labelid = "model", varid = "amp.avg", outlabel = "model"))

pdf(file.path(plot.dir, "summarize_livkid_modules.pdf"))
mm <- ggplot(fits.counts.by.amp, aes(x = 2 * amp.thres, y = n.genes, group = model, colour = as.factor(model))) + geom_line() +
  geom_line(size = 2) +
  theme_bw(20) +
  labs(colour = "# Rhythmic\nTissues") +
  theme(aspect.ratio=1,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab("Avg Amplitude of Rhythmic Tissues") + ylab("# Genes") + xlim(c(0.15, 6)) +
  scale_y_log10(breaks = c(1, 10, 100, 1000)) +
  geom_vline(xintercept = 2.8, linetype = "dotted") +
  scale_colour_brewer(palette = "Spectral")
print(mm)
print(mm + theme(legend.position = "none"))
dev.off()

# Single TFs underlie clock-independent liver-specific diurnal gen --------

dat.mean.wtko <- dat.wtko.collapsed %>%
  group_by(gene, tissue, geno) %>%
  summarise(exprs.mean = mean(exprs))

# jmod <- "Liver_SV129,Liver_BmalKO"
jmods <- c("Liver_SV129,Liver_BmalKO", "Liver_BmalKO", "Kidney_SV129", "Liver_SV129", "Liver_SV129,Kidney_SV129,Liver_BmalKO,Kidney_BmalKO")
for (jmod in jmods){
  print(paste("Finding regulators for:", jmod))
  genes.mod <- unique(as.character(subset(fits.long.filt, model == jmod)$gene))
  
  maraoutdir <- paste0("/home/yeung/projects/tissue-specificity/results/MARA.liver_kidney/promoters.", jmod, ".g=1001/atger_with_kidney.bugfixed")
  act.s <- LoadActivitiesLongKidneyLiver(maraoutdir, shorten.motif.name = TRUE)
  act.s.complex <- subset(act.s, tissue != "WFAT") %>%
    group_by(gene, tissue) %>%
    do(ProjectToFrequency2(., omega, add.tissue=TRUE))
  s.act <- SvdOnComplex(act.s.complex, value.var = "exprs.transformed")
  
  
  jmod.label <- gsub(pattern = ",", replacement = "-", jmod)
  pdf(file.path(plot.dir, paste0("clock_independent_", jmod.label, "_regulators.pdf")))
  
  s <- SvdOnComplex(subset(dat.freq, gene %in% genes.mod), value.var = "exprs.transformed")
  eigens <- GetEigens(s, period = 24, comp = comp, label.n = 20, eigenval = TRUE, adj.mag = TRUE, constant.amp = 4, peak.to.trough = TRUE, label.gene = c("Jun", "Egr1"))
  jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)
  
  eigens.act <- GetEigens(s.act, period = 24, comp = comp, adj.mag = TRUE, constant.amp = 4, label.n = 20, jtitle = "", peak.to.trough = TRUE)
  
  # plot genes and regulators
  print(eigens$u.plot)
  multiplot(eigens$u.plot, eigens$v.plot, layout = jlayout)
  print(eigens.act$u.plot)
  multiplot(eigens.act$u.plot, eigens.act$v.plot, cols = 2)
  # show mean exprs across tissues
  print(PlotMeanExprsOfModel(dat.mean.wtko, genes = genes.mod, jmodel = jmod, sorted = TRUE, avg.method = "mean"))
  
  # Print genes that match model
  jmotifs <- names(head(eigens.act$eigensamp[order(abs(eigens.act$eigensamp), decreasing = TRUE)], n = 20))
  genes.all <- unlist(sapply(jmotifs, GetGenesFromMotifs, tfs))
  genes.that.fit <- as.character(subset(fits.long.filt, gene %in% genes.all & model %in% jmod & method == jmeth)$gene)
  if (length(genes.that.fit) > 0){
    for (gene.hit in genes.that.fit){
      print(PlotGeneTissuesWTKO(subset(dat.wtko, gene == gene.hit), jtitle = gene.hit))
    }
  }
  dev.off()  
}


# Cooperative TFs underlie clock-dependent tissue-specific diurnal --------

# show enrichment of DHS peaks
wtmodulef <- "/home/yeung/projects/tissue-specificity/Robjs/liver_kidney_atger_nestle/tissue_specific_peaks/n.tiss.spec.df.out.lst.rand.1000.tissue.Liver.module.Liver_SV129random.flat.TRUE.Robj"
load(wtmodulef, v=T)
df.out.lst.merged.liverWT <- df.out.lst.merged
df.out.lst.merged.liverWT$gene.type[1:3] <- c("Liver_SV129", "Flat_SV129", "Flat.filt_SV129")

wtko.modulef <- "/home/yeung/projects/tissue-specificity/Robjs/liver_kidney_atger_nestle/tissue_specific_peaks/n.tiss.spec.df.out.lst.rand.1000.tissue.Liver.module.Liver_SV129,Liver_BmalKOrandom.flat.TRUE.Robj"
load(wtko.modulef, v=T)
df.out.lst.merged.liverWTKO <- df.out.lst.merged
df.out.lst.merged.liverWTKO$gene.type[1:3] <- c("Liver_SV129.Liver_BmalKO", "Flat_SV129BmalKO", "Flat.filt_SV129BmalKO")
df.out.lst.merged <- rbind(df.out.lst.merged.liverWT, df.out.lst.merged.liverWTKO)
df.out.lst.merged$xlabs <- make.names(df.out.lst.merged$gene.type, unique = TRUE)
df.out.lst.merged <- df.out.lst.merged[order(df.out.lst.merged$gene.type, decreasing = FALSE), ]

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
m.bar <- ggplot(df.out.lst.meanvar, aes(x = gene.type, y = mean.frac)) + geom_bar(stat = "identity", width = barwidth) + theme_bw() + 
  geom_errorbar(limits, width = barwidth / 2) + theme(aspect.ratio=1.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  
  xlab("") + ylab("Fraction of Genes with Liver-specific DHS peaks")



jmod <- "Liver_SV129,Liver_BmalKO"
load("Robjs/penalized_lda_robjs/mat.pmd.RORA_bHLH_SRF.40000.g1001.Robj", v=T)
pldarobj <- "/home/yeung/projects/tissue-specificity/Robjs/liver_kidney_atger_nestle/plda_robjs/2D.posterior.multigene.distfilt.morenonliv.bugfixed.liverWTKO.40000.cutoff.3.cutofflow0.method.g=1001.Robj"
load(pldarobj)

mat.fgbg.lab.lst.3 <- SetUpMatForLda(mat.fg, mat.bgnonliver, mat.bg, has.peaks = TRUE)
mat.fgbg.3 <- mat.fgbg.lab.lst.3$mat.fgbg; labels3 <- mat.fgbg.lab.lst.3$labels
colnames(mat.fgbg.3) <- sapply(colnames(mat.fgbg.3), RemoveP2Name)
colnames(mat.fgbg.3) <- sapply(colnames(mat.fgbg.3), function(cname){
  return(RemoveCommasBraces(cname))
}, USE.NAMES = FALSE)


rhyth.motifs <- sapply(GetTopMotifs("rhythmic"), RemoveP2Name, USE.NAMES = FALSE)
rhyth.motifs <- c(rhyth.motifs, c("SRF"))
# add DBP
# rhyth.motifs <- c(rhyth.motifs, c("DBP"))
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

jlambda <- 0.035  # liv only
out.cross.rhythtiss3 <- PenalizedLDA(mat.fgbg.cross.rhythtiss3, labels3, lambda = jlambda, K = 2, standardized = FALSE)

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
  geom_text_repel(data = dat.labs, aes(x = x, y = y, label = motif), size = 1.5) + 
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + theme_bw() + 
  theme(aspect.ratio = 0.33, legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  xlab("Motif loadings separating liver and kidney DHS peaks") + ylab("Motif loadings separating rhythmic and flat DHS peaks")


# MARA on subset of motifs

descrip <- "sep_liv_rhyth.RORA_bHLH_SRF.fixsignbugfixed"
K <- 3

outmain <- paste0("/home/yeung/projects/tissue-specificity/results/MARA.liver_kidney/mara_out.", descrip, ".K_", K)
inmain <- outmain
indir <- file.path(inmain, "atger_with_kidney.bugfixed")
# indir <- "/home/yeung/projects/tissue-specificity/results/MARA.liver_kidney/many_motifs_decorrelated.nofiltergenes.noabs/atger_with_kidney.bugfixed"
act.long <- LoadActivitiesLongKidneyLiver(indir, collapse.geno.tissue = TRUE, shorten.motif.name = FALSE)
# rename motifs based on the motifs with non-zero entries
act.long$pc <- sapply(as.character(act.long$gene), function(g) as.numeric(strsplit(g, "_")[[1]][[3]]), USE.NAMES = FALSE)
motif.decor <- sapply(act.long$pc, function(p) PcToMotif(mat.pmd$v, p), USE.NAMES = FALSE)
act.long$gene <- motif.decor

act.complex <- act.long %>%
  group_by(gene, tissue) %>%
  do(ProjectToFrequency2(., omega, add.tissue=TRUE))


s.act <- SvdOnComplex(act.complex, value.var = "exprs.transformed")

jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)
# jtitle <- gsub(pattern = "\\.", replacement = "\n", basename(indirmain))

max.labs <- 2
jtitle <- ""
eigens.act <- GetEigens(s.act, period = 24, comp = comp, adj.mag = TRUE, constant.amp = 5, label.n = K, jtitle = jtitle, peak.to.trough = TRUE)

pdf(file.path(plot.dir, "clock_dependent_liver_genes_cooperative_action.pdf"))
  print(m.bar)
  print(m)
  print(eigens.act$u.plot)
  multiplot(eigens.act$u.plot, eigens.act$v.plot, cols = 2)
dev.off()


# 4CSeq reveals tissue-specific promoter-enhancer contacts ----------------

projdir <- "/home/yeung/projects/4c_seq"
load(file.path(projdir, "Robjs/counts.long.delt.lk.LR.Robj"), v=T)
load(file.path(projdir, "Robjs/counts.hits.ignore5frags.normInf.with_assay_signal.Robj"), v=T)

source(file.path(projdir, "scripts/functions/GetBaitLocations.R"))
source(file.path(projdir, "scripts/functions/PlotFunctions.R"))
source(file.path(projdir, "scripts/functions/TrackHubFunctions.R"))
source(file.path(projdir, "scripts/functions/AssayFunctions.R"))
source(file.path(projdir, "scripts/functions/FindPeaks.R"))
source(file.path(projdir, "scripts/functions/MergeCounts.R"))

bait.locs <- GetBaitLocations()

zt <- "ZT20"
zt08 <- "ZT08"
jdist <- 250000

jbaits <- c("Mreg", "GYS2", "Pi3kap1", "Slc44a1", "Slc45a3short", "Slc45a3long")
baits.newnames <- c("Mreg", "Gys2", "Pik3ap1", "Slc44a1", "Slc45a3 (short)", "Slc45a3 (long)")
baits.genename <- c("Mreg", "Gys2", "Pik3ap1", "Slc44a1", "Slc45a3", "Slc45a3")
# jbaits <- c("Mreg", "Pi3kap1", "Slc44a1", "Slc45a3short", "Slc45a3long")
# baits.newnames <- c("Mreg", "Pik3ap1", "Slc44a1", "Slc45a3", "Slc45a3")
baits.hash <- hash(jbaits, baits.newnames)
baits.genename.hash <- hash(jbaits, baits.genename)

remove.bait <- c("Gys2", "Slc45a3 (long)")
# rename for volcano plots
counts.delt.lk.sub$bait <- factor(sapply(as.character(counts.delt.lk.sub$bait), function(b) baits.hash[[b]]), levels = unique(baits.newnames))

pdf(file.path(plot.dir, "4cseq_tissue_plots.nogys2.pdf"))

  for (jbait in jbaits){
    jbait.new <- baits.hash[[jbait]]
    jbait.genename.new <- baits.genename.hash[[jbait]]
    print(PlotGeneTissuesWTKO(subset(dat.wtko, gene == jbait.genename.new), jtitle = jbait.genename.new))
    
    jsub.sig <- subset(counts.long, genotype == "WT" & time == zt)
    jsub.sig08 <- subset(counts.long, genotype == "WT" & time == zt08)
    Signal <- PlotSignalLivVsKidLR(jbait, jsub.sig, pseudo.low = 500, jtitle = paste(jbait.new, zt), mindist = jdist, do.facet = FALSE)
    Signal08 <- PlotSignalLivVsKidLR(jbait, jsub.sig08, pseudo.low = 500, jtitle = paste(jbait.new, zt08), mindist = jdist, do.facet = FALSE, show.legend = FALSE)
    Zscores <- MergeCountsLivKid(jbait, counts.delt.lk, bait.locs, max.dist = jdist, show.plot = "Zscore", jtitle = "", jxlab = "", jshow.legend = FALSE, flip.y.axis = TRUE)
    Pvalues <- MergeCountsLivKid(jbait, counts.delt.lk, bait.locs, max.dist = jdist, show.plot = "Pvalue", jtitle = "", jshow.legend = TRUE, flip.y.axis = TRUE)
    PlotQuad(Signal, Signal08, Zscores, Pvalues, n.boxes = 2)
  }

  ggplot(subset(counts.delt.lk.sub, bait == "Mreg"), aes(x = A.delta / log10(2), y = -log10(pval.row))) +
    theme_bw() + 
    geom_point(alpha = 0.5) + facet_grid(time ~ bait) + xlim(-2.4, 2.4) + 
    xlab("Log fold change (Kidney - Liver)") + 
    ylab("-log10(p-value") + 
    theme(aspect.ratio = 1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  ggplot(subset(counts.delt.lk.sub, bait %in% baits.newnames & ! bait %in% remove.bait), aes(x = A.delta / log10(2), y = -log10(pval.row))) +
    theme_bw() + 
    geom_point(alpha = 0.5) + facet_grid(time ~ bait) + xlim(-2.4, 2.4) + 
    xlab("Log fold change (Kidney - Liver)") + 
    ylab("-log10(p-value") + 
    theme(aspect.ratio = 1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  ggplot(subset(counts.delt.lk.sub, bait %in% baits.newnames), aes(x = A.delta / log10(2), y = -log10(pval.row))) +
    theme_bw() + 
    geom_point(alpha = 0.5) + facet_grid(time ~ bait) + xlim(-2.4, 2.4) + 
    xlab("Log fold change (Kidney - Liver)") + 
    ylab("-log10(p-value") + 
    theme(aspect.ratio = 1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  

dev.off()