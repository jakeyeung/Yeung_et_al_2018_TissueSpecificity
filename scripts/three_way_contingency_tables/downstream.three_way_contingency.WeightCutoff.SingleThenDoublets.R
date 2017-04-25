# 2017-04-20
# Jake Yeung
# Use Weight cutoff and then run three way contingency table
# What makes FOX CUX ONECUT interesting? Or are others more interesting?

rm(list=ls())

library(dplyr)
library(hash)
library(ggplot2)
library(reshape2)

source("scripts/functions/PlotFunctions.R")
source("scripts/functions/PlotGeneAcrossTissues.R")

RunPoissonModel.test <- function(dat){
  mod1 <- glm(freq ~ model + motif1 * motif2, data=dat, family=poisson())
  chisq.pval <- 1 - pchisq(mod1$deviance, mod1$df.residual)
  return(data.frame(deviance = mod1$deviance, pval = chisq.pval))
}

RunPoissonModel.JI <- function(dat){
  mod1 <- glm(freq ~ model + motif1 + motif2 + motif1 * motif2, data=dat, family=poisson())
  chisq.pval <- 1 - pchisq(mod1$deviance, mod1$df.residual)
  return(data.frame(deviance = mod1$deviance, pval = chisq.pval))
}

RunPoissonModel.sing <- function(dat){
  mod1 <- glm(freq ~ model + motif1, data=dat, family=poisson())
  chisq.pval <- 1 - pchisq(mod1$deviance, mod1$df.residual)
  return(data.frame(deviance = mod1$deviance, pval = chisq.pval))
}

RunPoissonModel.cnd <- function(dat){
  mod1 <- glm(freq ~ model + motif1 + motif2 + model * motif1 + model * motif2, data=dat, family=poisson())
  return(mod1)
  chisq.pval <- 1 - pchisq(mod1$deviance, mod1$df.residual)
  return(data.frame(deviance = mod1$deviance, pval = chisq.pval))
}

RunPoissonModel.cnd2 <- function(dat){
  mod1 <- glm(freq ~ (motif1 + motif2) * model, data=dat, family=poisson())
  chisq.pval <- 1 - pchisq(mod1$deviance, mod1$df.residual)
  return(data.frame(deviance = mod1$deviance, pval = chisq.pval))
}

GetTargetGenesFromPair <- function(N.mat.all, pair, model = "rhyth", model.cname = "model"){
  motif1 <- strsplit(pair, ";")[[1]][[1]]
  motif2 <- strsplit(pair, ";")[[1]][[2]]
  cnames <- c("gene", "peak", "model", motif1, motif2)
  jsub <- N.mat.all[, cnames]
  jsub <- dplyr::filter_(jsub, paste(paste(motif1, "==",  "'atop'"), "&", paste(motif2, "==",  "'atop'"),"&", model.cname, "==", paste0("'", model, "'")))
  return(jsub)
}


# Load --------------------------------------------------------------------

load("Robjs/liver_kidney_atger_nestle/dat.long.liverkidneyWTKO.bugfixed.Robj", v=T); dat.wtko <- dat.long; rm(dat.long)


K <- 300
jweight <- 0.8
nb.models <- 3
nb.models <- 2
inf <- paste0("/home/yeung/projects/tissue-specificity/Robjs/three_way_cooccurence/three.way.cooccurrence.bugfixed.nmodels.", nb.models, ".K.", K, ".weight.", jweight, ".MergePeaks.FALSE.nullmodel.JI.withNmatallNmatfreqs.RemoveZeroCounts.Robj")

inf <- "/home/yeung/projects/tissue-specificity/Robjs/three_way_cooccurence/three.way.cooccurrence.bugfixed.nmodels.2.K.300.weight.0.8.MergePeaks.FALSE.nullmodel.CI.withNmatallNmatfreqs.RemoveZeroCounts.Robj"

load(inf, v=T)  # load also Nmat and Nmatfreqs

jmethod <- "g=1001"
load("Robjs/liver_kidney_atger_nestle/fits.long.multimethod.filtbest.staggeredtimepts.bugfixed.Robj", v=T)
fits.best.orig <- fits.long.filt
fits.best <- fits.long.filt; rm(fits.long.filt)
fits.best <- subset(fits.best, method == jmethod)
fits.best.filt <- subset(fits.best, weight >= jweight)
fits.best.sub <- subset(fits.best.filt, model == "Liver_SV129")


N.mat.freqs.all <- N.mat.freqs

if (length(unique(N.mat.freqs.all$model)) == 3){
  # ignore Kidney
  # N.mat.freqs <- subset(N.mat.freqs, model %in% c("flat", "rhyth"))
  
  # consider Kidney as flat!
  N.mat.freqs.all$model.merge <- sapply(N.mat.freqs.all$model, function(m) ifelse(m == "rhyth", "rhyth", "flat"))
  N.mat.freqs <- N.mat.freqs.all %>%
    group_by(model.merge, pair, motif1, motif2) %>%
    summarise(freq = sum(freq))
  N.mat.freqs <- dplyr::rename(N.mat.freqs, "model" = model.merge)
}



# Perform on singlets -----------------------------------------------------

N.mat.long <- melt(N.mat.all, id.vars = c("gene", "peak", "model"), variable.name = "motif")
N.mat.sum <- N.mat.long %>%
  group_by(model, motif) %>%
  filter(value == "atop") %>%
  filter(!duplicated(gene)) %>%
  summarise(n.genes.with.motif = length(which(value == "atop"))) %>%
  filter(model == "rhyth") %>%
  arrange(desc(n.genes.with.motif))  # ROR is most enriched motif!
N.mat.sum <- OrderDecreasing(N.mat.sum, "motif", "n.genes.with.motif")

N.mat.freqs.sing <- N.mat.long %>%
  group_by(model, motif, value) %>%
  summarise(freq = length(gene))
colnames(N.mat.freqs.sing) <- c("model", "motif", "motif1", "freq")


fits.single <- subset(N.mat.freqs.sing) %>%
  group_by(motif) %>%
  do(RunPoissonModel.sing(.)) %>%
  arrange(pval)

N.mat.OR.sing <- N.mat.freqs.sing %>%
  group_by(motif) %>%
  summarise(OR = (freq[3] / freq[4]) / (freq[1] / freq[2])) %>%
  arrange(desc(OR))

OR.sing <- hash(as.character(N.mat.OR.sing$motif), N.mat.OR.sing$OR)

fits.single$OR <- sapply(as.character(fits.single$motif), function(p) OR.sing[[p]])
fits.single$motif.label <- mapply(function(motif, p) ifelse(p < 0.05, motif, NA), as.character(fits.single$motif), fits.single$pval)

# Plot singlets -----------------------------------------------------------


# show enrichment of ROR
fits.single$motif.label <- gsub("RORA", "RORE", fits.single$motif.label)
m.single <- ggplot(fits.single, aes(x = log10(OR), y = -log10(pval), label = motif.label, size = -log10(pval))) + geom_point(size = 1) + geom_text_repel(nudge_y = 0.05) + 
  theme_bw(24) + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# show Amplitude and Phase of ROR target genes
load("Robjs/liver_kidney_atger_nestle/fits.bytiss.bugfixed.Robj", v=T); fits.bytiss <- subset(fits.bytiss, tissue == "Liver_SV129")
genes.liv.weightfilt <- as.character(fits.best.sub$gene)
load("Robjs/liver_kidney_atger_nestle/N.long.all_genes.3wtmodules.bugfixed.Robj", v=T); N.long.filt <- subset(N.long.livertwflat, motif == "RORA.p2"); rm(N.long.livertwflat)
peaks.liv <- as.character(subset(N.mat.all, model == "rhyth")$peak)  # 313 peaks
genes.liv <- unique(as.character(subset(N.mat.all, model == "rhyth")$gene))  # 139 genes
genes.liv.weightfilt <- as.character(fits.best.sub$gene)
N.sub <- subset(N.long.filt, peak %in% peaks.liv & gene %in% genes.liv & motif == "RORA.p2") %>%
  group_by(gene) %>%
  summarise(sitecount = sum(sitecount))
genes.liv.filt <- as.character(N.sub$gene)

top.n <- 50
ROR.counts.hash <- hash(as.character(N.sub$gene), N.sub$sitecount)
fits.bytiss.sub <- subset(fits.bytiss, gene %in% genes.liv.filt)
fits.bytiss.sub$ROR.sc <- sapply(as.character(fits.bytiss.sub$gene), function(g) ROR.counts.hash[[g]])
fits.bytiss.sub <- fits.bytiss.sub %>% arrange(desc(ROR.sc))

ror.hits <- as.character(fits.bytiss.sub[1:top.n, ]$gene)
fits.bytiss$is.hit <- sapply(as.character(fits.bytiss$gene), function(g) ifelse(g %in% ror.hits, TRUE, FALSE))
fits.bytiss$is.hit <- factor(fits.bytiss$is.hit, levels = c(TRUE, FALSE))

fits.bytiss.sum <- subset(fits.bytiss, gene %in% genes.liv.weightfilt) %>% 
  group_by(is.hit) %>% 
  summarise(n.genes = length(gene))

pdf(paste0("/home/yeung/projects/tissue-specificity/plots/liver_dhs_enrichment_plots/figure_5_singlets.pdf"), useDingbats = FALSE)
print(m.single)
ggplot(subset(fits.bytiss, gene %in% genes.liv.weightfilt), aes(x = is.hit, y = 2 * amp)) + geom_violin() + xlab("Has ROR Motifs") + ylab("Log2 Fold Change") + ggtitle("With RORE vs Liver_SV129 genes with BIC threshold") + 
  theme_bw(24) + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  annotate(geom = "text", x = 1, y = 4, label = paste("n=", fits.bytiss.sum$n.genes[[1]])) + 
  annotate(geom = "text", x = 2, y = 4, label = paste("n=", fits.bytiss.sum$n.genes[[2]]))

ggplot(subset(fits.bytiss, gene %in% genes.liv.weightfilt), aes(x = is.hit, y = 2 * amp)) + geom_boxplot() + xlab("Has ROR Motifs") + ylab("Log2 Fold Change") + ggtitle("With RORE vs Liver_SV129 genes with BIC threshold") + 
  theme_bw(24) + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  annotate(geom = "text", x = 1, y = 4, label = paste("n=", fits.bytiss.sum$n.genes[[1]])) + 
  annotate(geom = "text", x = 2, y = 4, label = paste("n=", fits.bytiss.sum$n.genes[[2]]))

library(plotrix)
jcex <- 2
circular_phase24H_histogram(subset(fits.bytiss, gene %in% genes.liv.weightfilt & is.hit == TRUE)$phase, jtitle = paste(""), cex.axis = jcex, cex.lab = jcex)
circular_phase24H_histogram(subset(fits.bytiss, gene %in% genes.liv.weightfilt & is.hit == FALSE)$phase, jtitle = paste(""), cex.axis = jcex, cex.lab = jcex)
dev.off()

# Plot doublets -----------------------------------------------------------

inf <- "/home/yeung/projects/tissue-specificity/Robjs/three_way_cooccurence/three.way.cooccurrence.bugfixed.nmodels.3.K.300.weight.0.MergePeaks.FALSE.nullmodel.JI.withNmatallNmatfreqs.RemoveZeroCounts.MergeOnecutCux.Robj"
inf <- "/home/yeung/projects/tissue-specificity/Robjs/three_way_cooccurence/three.way.cooccurrence.bugfixed.nmodels.3.K.300.weight.0.MergePeaks.FALSE.nullmodel.JI.MergeOnecutCux.TRUE.withNmatallNmatfreqs.RemoveZeroCounts.Robj"
inf <- "/home/yeung/projects/tissue-specificity/Robjs/three_way_cooccurence/three.way.cooccurrence.bugfixed.nmodels.3.K.300.weight.0.8.MergePeaks.FALSE.withNmatallNmatfreqs.RemoveZeroCounts2.Robj"
flatampmax <- 0.1
flatampmax <- 0.15
flatampmax <- 0.05
inf <- paste0("/home/yeung/projects/tissue-specificity/Robjs/three_way_cooccurence/three.way.cooccurrence.bugfixed.nmodels.3.K.300.weight.0.MergePeaks.FALSE.nullmodel.JI.flatampmax.", flatampmax, ".withNmatallNmatfreqs.RemoveZeroCounts.Robj")

load(inf, v=T)
N.mat.freqs.all <- N.mat.freqs

# consider Kidney as background
N.mat.freqs.all$model.merge <- sapply(N.mat.freqs.all$model, function(m) ifelse(m == "rhyth", "rhyth", "flat"))
N.mat.freqs <- N.mat.freqs.all %>%
  group_by(model.merge, pair, motif1, motif2) %>%
  summarise(freq = sum(freq))
N.mat.freqs <- dplyr::rename(N.mat.freqs, "model" = model.merge)


# calculate some enrichment
ratios <- N.mat.freqs %>%
  group_by(pair, model) %>%
  summarise(odds.ratio = (freq[1] * freq[4]) / (freq[2] * freq[3]))
ratios.all <- N.mat.freqs.all %>%
  group_by(pair, model) %>%
  summarise(odds.ratio = (freq[1] * freq[4]) / (freq[2] * freq[3]))
ratios.TopVBottom <- N.mat.freqs %>%
  group_by(pair, model) %>%
  summarise(TvB = (freq[1] / freq[4]))
ratios.ofTvB <- ratios.TopVBottom %>%
  group_by(pair) %>%
  summarise(TvB.TvB = TvB[2] / TvB[1])
ratios.TopVBottom.rhyth <- subset(ratios.TopVBottom, model == "rhyth")
ratios.rhyth <- subset(ratios, model == "rhyth")
ratios.ofOR <- ratios %>%
  group_by(pair) %>%
  summarise(OR.OR = odds.ratio[2] / odds.ratio[1])

# add to fits
enrich <- hash(ratios.rhyth$pair, ratios.rhyth$odds.ratio)
enrich.OR.OR <- hash(ratios.ofOR$pair, ratios.ofOR$OR.OR)
enrich.TvB <- hash(ratios.TopVBottom.rhyth$pair, ratios.TopVBottom.rhyth$TvB)
enrich.TvB.TvB <- hash(ratios.ofTvB$pair, ratios.ofTvB$TvB.TvB)
enrich.NTop <- hash(subset(N.mat.freqs, model == "rhyth" & motif1 == "atop" & motif2 == "atop")$pair, subset(N.mat.freqs, model == "rhyth" & motif1 == "atop" & motif2 == "atop")$freq)

fits$OR <- sapply(fits$pair, function(p) enrich[[p]])
fits$OR.OR <- sapply(fits$pair, function(p) enrich.OR.OR[[p]])
fits$TvB <- sapply(fits$pair, function(p) enrich.TvB[[p]])
fits$TvB.TvB <- sapply(fits$pair, function(p) enrich.TvB.TvB[[p]])
fits$NTop <- sapply(fits$pair, function(p) enrich.NTop[[p]])

# pdf("/home/yeung/projects/tissue-specificity/plots/liver_dhs_enrichment_plots/figure_5_doublets.pdf", useDingbats = FALSE)
pdf(paste0("/home/yeung/projects/tissue-specificity/plots/liver_dhs_enrichment_plots/figure_5_doublets.flatampmax.", flatampmax, ".pdf"), useDingbats = FALSE)

jmotif <- "RORE"
fits$pair <- gsub("RORA", "RORE", fits$pair)

ggplot(subset(fits, grepl("RORE", fits$pair)), aes(x = log10(OR.OR), y = -log10(pval), label = pair)) + geom_point() + geom_text(size = 3) + geom_vline(xintercept = 0) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("Log10 Odds Ratio Compared to Background")

jfits <- subset(fits, grepl("RORE", fits$pair) & log10(OR.OR) > 0)
jfits$label <- mapply(function(pair, pval) ifelse(-log10(pval) > 6.5, pair, NA), jfits$pair, jfits$pval)
ggplot(jfits, aes(x = log10(OR.OR), y = -log10(pval), label = label)) + geom_point() + geom_text_repel(size = 6.5) + geom_vline(xintercept = 0) + 
  theme_bw(24) + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
  xlab("Log10 Odds Ratio Compared to Background")

jfits <- subset(fits, grepl("RORE", fits$pair) & log10(OR.OR) > 0)
jfits$label <- mapply(function(pair, pval) ifelse(-log10(pval) > 6.5, pair, NA), jfits$pair, jfits$pval)
ggplot(jfits, aes(x = OR.OR, y = -log10(pval), label = label)) + geom_point() + geom_text_repel(size = 6.5) + geom_vline(xintercept = 1) + 
  theme_bw(24) + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
  xlab("Odds Ratio Compared to Background")
dev.off()


# Plot enrichment of fraction of genes ------------------------------------

wtmodulef <- "/home/yeung/projects/tissue-specificity/Robjs/liver_kidney_atger_nestle/tissue_specific_peaks/n.tiss.spec.df.out.lst.rand.1000.tissue.Liver.module.Liver_SV129.random.flat.TRUE.BICweightCutoff.0.8.Robj"
load(wtmodulef, v=T)
df.out.lst.merged.liverWT <- df.out.lst.merged
df.out.lst.merged.liverWT$gene.type[1:3] <- c("Liver_SV129", "Flat_SV129", "Flat.filt_SV129")

wtko.modulef <- "/home/yeung/projects/tissue-specificity/Robjs/liver_kidney_atger_nestle/tissue_specific_peaks/n.tiss.spec.df.out.lst.rand.1000.tissue.Liver.module.Liver_SV129,Liver_BmalKO.random.flat.TRUE.BICweightCutoff.0.8.Robj"
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

# create coloring
df.out.lst.meanvar$model <- factor(c("Liver", "Liver", "Flat"), levels = c("Liver", "Flat"))
df.out.lst.meanvar$geno <- factor(c("WT", "Bmal1 KO", "WT"), levels = c("WT", "Bmal1 KO"))
df.out.lst.meanvar$gene.type <- factor(c("Clock-Driven", "System-Driven", "Nonrhythmic"), levels = c("Clock-Driven", "System-Driven", "Nonrhythmic"))

jcols <- c(scales::hue_pal()(2)[[1]], "gray")  # red (matching if two comparisons) and gray

# pdf("/home/yeung/projects/tissue-specificity/plots/liver_dhs_enrichment_plots/figure_5_fraction_liver_dhs_plots_colored.pdf", useDingbats = FALSE)
pdf(paste0("/home/yeung/projects/tissue-specificity/plots/liver_dhs_enrichment_plots/figure_5_liver_dhs_plots_colored.flatampmax.pdf"), useDingbats = FALSE)
barwidth <- 0.8
limits <- aes(ymax = mean.frac + sd.frac, ymin=mean.frac - sd.frac)
m.bar <- ggplot(df.out.lst.meanvar, aes(x = gene.type, y = mean.frac, fill = model, linetype = geno)) + geom_bar(stat = "identity", width = barwidth, size = 1, colour = "black") + 
  theme_bw(24) + 
  geom_errorbar(limits, width = barwidth / 2) + 
  # theme(aspect.ratio=1.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  
  theme(axis.text.x = element_text(angle = 45, hjust = 1), aspect.ratio=1.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  scale_fill_manual(values = jcols) +
  ylab("Fraction of Genes") + xlab("")
print(m.bar)
dev.off()


# Activity of Liver TFs across DHS ----------------------------------------

source("scripts/functions/LoadActivitiesLong.R")
source("scripts/functions/PlotActivitiesFunctions.R")
source("scripts/functions/SortByTissue.R")

indir <- "/home/yeung/projects/tissue-specificity/results/MARA/MARA_DHS/MARA_dhs_gene_region_cutoff_filter/"
zscores.path <- "/home/yeung/projects/tissue-specificity/results/MARA/MARA_DHS/MARA_dhs_gene_region_cutoff_filter/"
act.dhs.long <- LoadActivitiesLongDhs(indir, act.file = "Activities", se.file = "StandardError")
zscores <- read.table(zscores.path, header=FALSE, sep="\t", col.names = c("motif", "zscore"))
zscores <- zscores[order(zscores$zscore, decreasing = TRUE), ]

act.sum <- act.dhs.long %>%
  group_by(gene) %>%
  summarise(exprs.sum = sum(exprs))

hits <- c("ONECUT1.2.p2", "NKX2.1.4.p2", "HNF4A_NR2F1.2.p2", "FOXA2.p3", "RFX1..5_RFXANK_RFXAP.p2", "CUX2.p2", "HNF1A.p2", "TEAD1.p2", "TFAP2.A.C..p2")
hits.rnaseq <- c("Onecut1", "Nkx2-1", "Hnf4a", "Foxa2", "Rfxank", "Onecut1", "Hnf1a", "Tead1", "Tfap2a")

# Load RNA-Seq from ENCODE data to plot side by side
load("Robjs/dat.encode.Robj")

dhsplots <- lapply(hits, function(jgene) {
  dat.sub <- SortByTissue(subset(act.dhs.long, gene == jgene), by.var = "exprs")
  return(PlotActivitiesWithSE.dhs(SortByTissue(subset(act.dhs.long, gene == jgene), by.var = "exprs")))
})

encodeplots <- lapply(hits.rnaseq, function(jgene) {
  dat.sub <- SortByTissue(subset(dat.encode, gene == jgene & !tissue %in% c("Colon", "Duodenum", "Small Intestine", "Stomach", "Placenta")), by.var = "tpm")
  m <- ggplot(dat.sub, aes(x = tissue, y = tpm)) + geom_bar(stat = "identity") + ggtitle(jgene) + xlab("") + ylab("TPM")
  return(m)
})

for (i in 1:length(hits)){
  dhsgene <- hits[i]
  encodegene <- hits.rnaseq[i]
  dhsplot <- PlotActivitiesWithSE.dhs(subset(act.dhs.long, gene == dhsgene))
  encodeplot <- PlotEncodeRnaseq(subset(dat.encode, gene == encodegene & !tissue %in% c("Colon", "Duodenum", "Small Intestine", "Stomach", "Placenta")))
  
  multiplot(dhsplot, encodeplot, cols = 2)
}


# Are ROR hits also rhythmic in lung? -------------------------------------

load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.15.phase_sd_maxdiff_avg.Robj", v=T)
load("Robjs/dat.long.fixed_rik_genes.Robj", v=T)  # hogenesch

foxa2.rora.hits <- GetTargetGenesFromPair(N.mat.all, pair = "FOXA2;RORA")

fits.livlung <- subset(fits.best, gene %in% ror.hits & grepl("Lung", model))

# pdf("/home/yeung/projects/tissue-specificity/plots/liver_dhs_enrichment_plots/foxa2_ror_hits_all_tissues.pdf", useDingbats = FALSE)
pdf(paste0("/home/yeung/projects/tissue-specificity/plots/liver_dhs_enrichment_plots/foxa2_rra-hits_all_tissues.flatampmax.", flatampmax, ".pdf"), useDingbats = FALSE)
for (g in as.character(foxa2.rora.hits$gene)){
  print(g)
  if (nrow(subset(dat.long, gene == g)) == 0){
    print(paste('Skipping', g))
    next
  }
  print(PlotGeneAcrossTissues(subset(dat.long, gene == g)))
}
dev.off()



# What are the flat hits? -------------------------------------------------

GetAmpOfFgBg <- function(fits.bytiss, N.mat.all, jpair){
  fg.hits <- GetTargetGenesFromPair(N.mat.all, pair = jpair, model = "rhyth", model.cname = "model")
  bg.hits <- GetTargetGenesFromPair(N.mat.all, pair = jpair, model = "flat", model.cname = "model")
  genes.fgbg <- c(as.character(fg.hits$gene), as.character(bg.hits$gene))
  jhash <- hash(genes.fgbg, c(rep(TRUE, length(fg.hits$gene)), rep(FALSE, length(bg.hits$gene))))
  jsub <- subset(fits.bytiss, gene %in% genes.fgbg)
  jsub$is.fg <- sapply(as.character(jsub$gene), function(g){
    is.fg <- jhash[[g]]
    if (is.null(is.fg)){
      warning(paste("Gene not in hash table:", g))
    }
    return(is.fg)
  })
  return(jsub)
}

jpair <- "ONECUT1.2;RORA"
jpair <- "CUX2;RORA"
jpair <- "FOXA2;RORA"
jpair <- "RORA;TFAP2.A.C."
jpair <- "RORA;HNF1A"
jsub <- GetAmpOfFgBg(fits.bytiss, N.mat.all, jpair = jpair)
ggplot(jsub, aes(x = amp, fill = is.fg)) + geom_density(alpha = 0.5) + ggtitle(jpair)
subset(jsub, amp < 0.15)

foxa2.rora.bg <- GetTargetGenesFromPair(N.mat.all, pair = "FOXA2;RORA", model = "flat", model.cname = "model")
onecut1.rora.bg <- GetTargetGenesFromPair(N.mat.all, pair = "ONECUT1.2;RORA", model = "flat", model.cname = "model")
cux2.rora.bg <- GetTargetGenesFromPair(N.mat.all, pair = "CUX2;RORA", model = "flat", model.cname = "model")

# Are flat hits "rhythmic"?
ggplot(subset(fits.bytiss, gene %in% foxa2.rora.bg$gene), aes(x = amp)) + geom_density()
ggplot(subset(fits.bytiss, gene %in% onecut1.rora.bg$gene), aes(x = amp)) + geom_density()
ggplot(subset(fits.bytiss, gene %in% cux2.rora.bg$gene), aes(x = amp)) + geom_density()
ggplot(subset(fits.bytiss, gene %in% foxa2.rora.hits$gene), aes(x = amp)) + geom_density()

genes.foxa2rora <- c(as.character(foxa2.rora.hits$gene), as.character(foxa2.rora.bg$gene))
foxa2rora.hash <- hash(genes.foxa2rora, c(rep(TRUE, length(foxa2.rora.hits$gene)), rep(FALSE, length(foxa2.rora.bg$gene))))
jsub <- subset(fits.bytiss, gene %in% genes.foxa2rora)
jsub$is.fg <- sapply(as.character(jsub$gene), function(g) foxa2rora.hash[[g]])
ggplot(jsub, aes(x = amp, fill = is.fg)) + geom_density(alpha = 0.5) + ggtitle("FOXA2;RORA")



# fits.bytiss$fg <- sapply(as.character(fits.bytiss$gene), function(g) ifelse(g %in% foxa2.rora.hits$gene, TRUE, FALSE))
# fits.bytiss$bg <- sapply(as.character(fits.bytiss$gene), function(g) ifelse(g %in% foxa2.rora.bg$gene, TRUE, FALSE))


# PlotGeneAcrossTissues(subset(dat.long, gene == "Insig2"))
# PlotGeneAcrossTissues(subset(dat.long, gene == "Slc4a4"))
# PlotGeneAcrossTissues(subset(dat.long, gene == "Abcg5"))
# PlotGeneAcrossTissues(subset(dat.long, gene == "Adrb3"))
# PlotGeneAcrossTissues(subset(dat.long, gene == "Arhgap18"))
# PlotGeneAcrossTissues(subset(dat.long, gene == ""))
