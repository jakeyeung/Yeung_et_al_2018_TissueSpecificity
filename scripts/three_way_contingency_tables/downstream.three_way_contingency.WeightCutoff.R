# 2017-04-20
# Jake Yeung
# Use Weight cutoff and then run three way contingency table
# What makes FOX CUX ONECUT interesting? Or are others more interesting?

rm(list=ls())

library(dplyr)
library(hash)
library(ggplot)

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

K <- 200
weight.cutoff <- 0.8
# inf <- "/home/yeung/projects/tissue-specificity/Robjs/three_way_cooccurence/three.way.cooccurrence.bugfixed.nmodels.2.K.300.withNmatallNmatfreqs.RemoveZeroCounts.Robj"
inf <- paste0("/home/yeung/projects/tissue-specificity/Robjs/three_way_cooccurence/three.way.cooccurrence.bugfixed.nmodels.2.K.", K, ".weight.", weight.cutoff, ".withNmatallNmatfreqs.RemoveZeroCounts.Robj")

inf <- "/home/yeung/projects/tissue-specificity/Robjs/three_way_cooccurence/three.way.cooccurrence.bugfixed.nmodels.2.K.200.weight.0.8.MergePeaks.TRUE.withNmatallNmatfreqs.RemoveZeroCounts.Robj"
inf <- "/home/yeung/projects/tissue-specificity/Robjs/three_way_cooccurence/three.way.cooccurrence.bugfixed.nmodels.2.K.300.weight.0.8.MergePeaks.FALSE.withNmatallNmatfreqs.RemoveZeroCounts.LiverPeaksvsBoth.Robj"
inf <- "/home/yeung/projects/tissue-specificity/Robjs/three_way_cooccurence/three.way.cooccurrence.bugfixed.nmodels.2.K.300.weight.0.8.MergePeaks.FALSE.withNmatallNmatfreqs.RemoveZeroCounts.LiverPeaksvsBoth.Robj"
inf <- "/home/yeung/projects/tissue-specificity/Robjs/three_way_cooccurence/three.way.cooccurrence.bugfixed.nmodels.2.K.300.weight.0.8.MergePeaks.FALSE.withNmatallNmatfreqs.RemoveZeroCounts2.Robj"

inf <- "/home/yeung/projects/tissue-specificity/Robjs/three_way_cooccurence/three.way.cooccurrence.bugfixed.nmodels.2.K.300.weight.0.8.MergePeaks.FALSE.nullmodel.CI.withNmatallNmatfreqs.RemoveZeroCounts.Robj"


inf <- "/home/yeung/projects/tissue-specificity/Robjs/three_way_cooccurence/three.way.cooccurrence.bugfixed.nmodels.2.K.300.weight.0.8.MergePeaks.FALSE.nullmodel.CI.withNmatallNmatfreqs.RemoveZeroCounts.Robj"


inf <- "/home/yeung/projects/tissue-specificity/Robjs/three_way_cooccurence/three.way.cooccurrence.bugfixed.nmodels.2.K.300.weight.0.8.MergePeaks.FALSE.withNmatallNmatfreqs.RemoveZeroCounts2.Robj"

inf <- "/home/yeung/projects/tissue-specificity/Robjs/three_way_cooccurence/three.way.cooccurrence.bugfixed.nmodels.3.K.300.weight.0.8.MergePeaks.FALSE.nullmodel.CI.withNmatallNmatfreqs.RemoveZeroCounts.Robj"

inf <- "/home/yeung/projects/tissue-specificity/Robjs/three_way_cooccurence/three.way.cooccurrence.bugfixed.nmodels.3.K.300.weight.0.8.MergePeaks.FALSE.withNmatallNmatfreqs.RemoveZeroCounts2.Robj"


K <- 300
weight <- 0.8
nb.models <- 3
inf <- paste0("/home/yeung/projects/tissue-specificity/Robjs/three_way_cooccurence/three.way.cooccurrence.bugfixed.nmodels.", nb.models, ".K.", K, ".weight.", weight, ".MergePeaks.FALSE.nullmodel.JI.withNmatallNmatfreqs.RemoveZeroCounts.Robj")

load(inf, v=T)  # load also Nmat and Nmatfreqs
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


# Combine CUX2 and ONECUT1 together? --------------------------------------



# Analyze -----------------------------------------------------------------


# Plot top hits by P-value AS WELL AS enrichment
# enrichment: pairs in rhyth over pairs in flat

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

ggplot(subset(fits, grepl("RORA", fits$pair)), aes(x = log10(OR), y = -log10(pval), label = pair)) + geom_point() + geom_text() + geom_vline(xintercept = 0) 
ggplot(subset(fits), aes(x = log10(OR), y = -log10(pval), label = pair)) + geom_point() + geom_text() + geom_vline(xintercept = 0) 
ggplot(subset(fits, grepl("RORA", fits$pair)), aes(x = log10(OR.OR), y = -log10(pval), label = pair)) + geom_point() + geom_text() + geom_vline(xintercept = 0) 
ggplot(subset(fits), aes(x = log10(OR.OR), y = -log10(pval), label = pair)) + geom_point() + geom_text() + geom_vline(xintercept = 0) 

ggplot(subset(fits, grepl("RORA", fits$pair) & pval < 1e-4), aes(x = log10(OR.OR), y = log10(TvB.TvB), label = pair)) + geom_point() + geom_text() + geom_vline(xintercept = 0) 
ggplot(subset(fits, grepl("RORA", fits$pair) & pval < 1e-4), aes(x = log10(OR.OR), y = NTop, label = pair)) + geom_point() + geom_text() + geom_vline(xintercept = 0) 

# Top single hits ---------------------------------------------------------

# Count top single hits for all motifs
N.mat.long <- melt(N.mat.all, id.vars = c("gene", "peak", "model"), variable.name = "motif")
N.mat.sum <- N.mat.long %>%
  group_by(model, motif) %>%
  filter(value == "atop") %>%
  filter(!duplicated(gene)) %>% 
  summarise(n.genes.with.motif = length(which(value == "atop"))) %>%
  filter(model == "rhyth") %>%
  arrange(desc(n.genes.with.motif))  # ROR is most enriched motif!
N.mat.sum <- OrderDecreasing(N.mat.sum, "motif", "n.genes.with.motif")

# 50 genes with ROR motif (313 peaks), that is amazing because you're competing with Kidney and Flat model
# Flat and Kidney has 1666 and 483 peaks respectively


ggplot(N.mat.sum[1:15, ], aes(x = motif, y = n.genes.with.motif)) + geom_bar(stat = "identity")

ggplot(subset(fits, grepl("RORA", fits$pair)), aes(x = log10(OR.OR), y = -log10(pval), label = pair)) + geom_point() + geom_text() + geom_vline(xintercept = 0) 
ggplot(subset(fits, grepl("FOXA2", fits$pair)), aes(x = log10(OR.OR), y = -log10(pval), label = pair)) + geom_point() + geom_text() + geom_vline(xintercept = 0) 
ggplot(subset(fits, grepl("ONECUT1", fits$pair)), aes(x = log10(OR.OR), y = -log10(pval), label = pair)) + geom_point() + geom_text() + geom_vline(xintercept = 0) 
ggplot(subset(fits, grepl("CUX2", fits$pair)), aes(x = log10(OR.OR), y = -log10(pval), label = pair)) + geom_point() + geom_text() + geom_vline(xintercept = 0) 

ggplot(subset(fits), aes(x = log10(OR.OR), y = -log10(pval), label = pair)) + geom_point() + geom_text() + geom_vline(xintercept = 0) 
ggplot(subset(fits, abs(log10(OR.OR)) > 0.3), aes(x = log10(OR.OR), y = -log10(pval), label = pair)) + geom_point() + geom_text() + geom_vline(xintercept = 0) 

ggplot(subset(fits, abs(log10(TvB.TvB)) > 0.3), aes(x = log10(TvB.TvB), y = -log10(pval), label = pair)) + geom_point() + geom_text() + geom_vline(xintercept = 0) 

ggplot(subset(fits, !grepl("ELF1.2", fits$pair)), aes(x = log10(OR.OR), y = -log10(pval), label = pair)) + geom_point() + geom_text() + geom_vline(xintercept = 0) 
ggplot(subset(fits, grepl("ELF1.2", fits$pair)), aes(x = log10(OR.OR), y = -log10(pval), label = pair)) + geom_point() + geom_text() + geom_vline(xintercept = 0) 
ggplot(subset(fits, grepl("HLF", fits$pair)), aes(x = log10(OR.OR), y = -log10(pval), label = pair)) + geom_point() + geom_text() + geom_vline(xintercept = 0) 

ggplot(subset(fits, grepl("RORA", fits$pair)), aes(x = TvB, y = -log10(pval), label = pair)) + geom_point() + geom_text()
ggplot(subset(fits), aes(x = TvB, y = -log10(pval), label = pair)) + geom_point() + geom_text()

ggplot(subset(fits, grepl("RORA", fits$pair)), aes(x = log10(TvB.TvB), y = -log10(pval), label = pair)) + geom_point() + geom_text()



# Plot primetime ----------------------------------------------------------

#jmotifs <- c("FOXA2", "ONECUT1.2", "CUX2")
jmotifs <- c("FOXA2")
jmotifs <- c("ONECUT1.2")
jmotifs <- c("RORA")
grepstr <- paste(jmotifs, collapse = "|")
(fits.sub <- subset(fits, grepl(grepstr, pair)) %>% arrange(pval))

fits.sub$has.ROR <- sapply(fits.sub$pair, function(p) grepl("RORA", p))
fits.sub$has.LivTF <- sapply(fits.sub$pair, function(p) grepl("FOXA2|ONECUT|CUX2", p))

fits.sub$pair <- factor(fits.sub$pair, levels = fits.sub$pair)

barcols <- c("grey50", "black")
m.loglin <- ggplot(subset(fits.sub, pval < 0.05 & OR.OR > 1)[1:20, ], aes(x = pair, y = -log10(pval), fill = has.LivTF)) + geom_bar(stat = "identity") +  
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle("Top motif pairs between\nFOXA2, ONECUT, or CUX2") + 
  scale_fill_manual(values = barcols) + xlab("")
print(m.loglin)

m.loglin2 <- ggplot(subset(fits.sub, pval < 0.05), aes(x = log10(OR.OR), y = -log10(pval), colour = has.ROR, label = pair)) + geom_point() +  geom_text() + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle("Top motif pairs between\nFOXA2, ONECUT, or CUX2") + 
  scale_fill_manual(values = barcols) + xlab("")
print(m.loglin2)



# Plot output: Amp Phase and Enrichment -----------------------------------



load("Robjs/liver_kidney_atger_nestle/N.long.all_genes.3wtmodules.bugfixed.Robj", v=T); N.long.filt <- N.long.livertwflat; rm(N.long.livertwflat)

weight.cutoff <- 0.8
jmethod <- "g=1001"
load("Robjs/liver_kidney_atger_nestle/fits.long.multimethod.filtbest.staggeredtimepts.bugfixed.Robj", v=T)
fits.best.orig <- fits.long.filt
fits.best <- fits.long.filt; rm(fits.long.filt)
fits.best <- subset(fits.best, method == jmethod)
fits.best.filt <- subset(fits.best, weight >= weight.cutoff)
fits.best.sub <- subset(fits.best.filt, model == "Liver_SV129")

load("Robjs/liver_kidney_atger_nestle/fits.bytiss.bugfixed.Robj", v=T)
fits.bytiss <- subset(fits.bytiss, tissue == "Liver_SV129")

genes.liv.all <- as.character(subset(fits.best.filt, model == "Liver_SV129")$gene)
genes.liv.weightfilt <- as.character(fits.best.sub$gene)
genes.liv <- unique(as.character(subset(N.mat.all, model == "rhyth")$gene))  # 139 genes
peaks.liv <- as.character(subset(N.mat.all, model == "rhyth")$peak)  # 313 peaks

# take top N genes with ROR elements, plot amplitude and phase
top.n <- 50
N.sub <- subset(N.long.filt, peak %in% peaks.liv & gene %in% genes.liv & motif == "RORA.p2") %>%
  group_by(gene) %>%
  summarise(sitecount = sum(sitecount))

N.sub.tissTF <- subset(N.long.filt, peak %in% peaks.liv & gene %in% genes.liv & motif %in% c("FOXA2.p3", "ONECUT1.2.p2", "HNF1A.p2", "CUX2.p2", "HNF4A_NR2F1,2.p2")) %>%
  group_by(gene) %>%
  summarise(sitecount = sum(sitecount),
            motif = "Liver") %>%
  group_by(gene)
N.sub.ROR <- subset(N.long.filt, peak %in% peaks.liv & gene %in% genes.liv & motif %in% c("RORA.p2")) %>%
  group_by(gene) %>%
  summarise(sitecount = sum(sitecount),
            motif = "RORA") %>%
  group_by(gene)

N.sub.LivROR <- bind_rows(N.sub.tissTF, N.sub.ROR) %>%
  group_by(gene) %>%
  filter(length(gene) == 2) %>%
  summarise(sitecount = prod(sitecount))


genes.liv.filt <- as.character(N.sub$gene)

# ROR
ROR.counts.hash <- hash(as.character(N.sub$gene), N.sub$sitecount)
fits.bytiss.sub <- subset(fits.bytiss, gene %in% genes.liv.filt)
fits.bytiss.sub$ROR.sc <- sapply(as.character(fits.bytiss.sub$gene), function(g) ROR.counts.hash[[g]])
fits.bytiss.sub <- fits.bytiss.sub %>% arrange(desc(ROR.sc))

# take top 50
ror.hits <- as.character(fits.bytiss.sub[1:top.n, ]$gene)
fits.bytiss$is.hit <- sapply(as.character(fits.bytiss$gene), function(g) ifelse(g %in% ror.hits, TRUE, FALSE))


# Take RORA + TF motifs
LivROR.counts.hash <- hash(as.character(N.sub.LivROR$gene), N.sub.LivROR$sitecount)
fits.bytiss.sub$LivROR.sc <- sapply(as.character(fits.bytiss.sub$gene), function(g){
  sc <- LivROR.counts.hash[[g]]
  return(ifelse(is.null(sc), NA, sc))
})
fits.bytiss.sub <- fits.bytiss.sub %>% arrange(desc(LivROR.sc))
livror.hits <- as.character(fits.bytiss.sub[1:top.n, ]$gene)
fits.bytiss$is.hit.livror <- sapply(as.character(fits.bytiss$gene), function(g) ifelse(g %in% livror.hits, TRUE, FALSE))

fits.bytiss$is.hit.livror <- factor(as.character(fits.bytiss$is.hit.livror), levels = c("TRUE", "FALSE"))
fits.bytiss$is.hit <- factor(as.character(fits.bytiss$is.hit.livror), levels = c("TRUE", "FALSE"))


# plot RORA enrichment with partners
pdf(paste0("/home/yeung/projects/tissue-specificity/plots/liver_dhs_enrichment_plots/figure_5_liver_dhs_plots.K.", K, ".nbmodels.", nb.models, ".weight.", weight, ".pdf"), useDingbats = FALSE)

ggplot(subset(fits, grepl("RORA", fits$pair)), aes(x = log10(OR.OR), y = -log10(pval), label = pair)) + geom_point() + geom_text(size = 3) + geom_vline(xintercept = 0) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("Log10 Odds Ratio Compared to Background")

jfits <- subset(fits, grepl("RORA", fits$pair) & log10(OR.OR) > 0)
jfits$label <- mapply(function(pair, pval) ifelse(-log10(pval) > 6.5, pair, NA), jfits$pair, jfits$pval)
ggplot(jfits, aes(x = log10(OR.OR), y = -log10(pval), label = label)) + geom_point() + geom_text_repel(size = 6.5) + geom_vline(xintercept = 0) + 
  theme_bw(24) + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
  xlab("Log10 Odds Ratio Compared to Background")

ggplot(jfits, aes(x = OR.OR, y = -log10(pval), label = label)) + geom_point() + geom_text_repel(size = 6.5) + geom_vline(xintercept = 1) + 
  theme_bw(24) + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
  xlab("Odds Ratio Compared to Background")

# amp

fits.bytiss.sum <- subset(fits.bytiss, gene %in% genes.liv.weightfilt) %>% 
  group_by(is.hit) %>% 
  summarise(n.genes = length(gene))

ggplot(subset(fits.bytiss, gene %in% genes.liv.weightfilt), aes(x = is.hit, y = 2 * amp)) + geom_violin() + xlab("Has ROR Motifs") + ylab("Log2 Fold Change") + ggtitle("With RORE vs Liver_SV129 genes with BIC threshold") + 
  theme_bw(24) + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  annotate(geom = "text", x = 1, y = 4, label = paste("n=", fits.bytiss.sum$n.genes[[1]])) + 
  annotate(geom = "text", x = 2, y = 4, label = paste("n=", fits.bytiss.sum$n.genes[[2]]))
  
ggplot(subset(fits.bytiss, gene %in% genes.liv.weightfilt), aes(x = is.hit.livror, y = 2 * amp)) + geom_violin() + xlab("Has ROR Motifs") + ylab("Log2 Fold Change") + ggtitle("With RORE vs Liver_SV129 genes with BIC threshold") + 
  theme_bw(24) + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  annotate(geom = "text", x = 1, y = 4, label = paste("n=", fits.bytiss.sum$n.genes[[1]])) + 
  annotate(geom = "text", x = 2, y = 4, label = paste("n=", fits.bytiss.sum$n.genes[[2]]))

# phase
library(plotrix)
jcex <- 2
circular_phase24H_histogram(subset(fits.bytiss, gene %in% genes.liv.weightfilt & is.hit == TRUE)$phase, jtitle = paste(""), cex.axis = jcex, cex.lab = jcex)
circular_phase24H_histogram(subset(fits.bytiss, gene %in% genes.liv.weightfilt & is.hit == FALSE)$phase, jtitle = paste(""), cex.axis = jcex, cex.lab = jcex)

circular_phase24H_histogram(subset(fits.bytiss, gene %in% genes.liv.weightfilt & is.hit.livror == TRUE)$phase, jtitle = paste(""), cex.axis = jcex, cex.lab = jcex)
circular_phase24H_histogram(subset(fits.bytiss, gene %in% genes.liv.weightfilt & is.hit.livror == FALSE)$phase, jtitle = paste(""), cex.axis = jcex, cex.lab = jcex)


dev.off()


# Plot Output Fraction of Genes with DHS: Weight Filtered -----------------


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


pdf("/home/yeung/projects/tissue-specificity/plots/liver_dhs_enrichment_plots/figure_5_fraction_liver_dhs_plots.pdf")
barwidth <- 0.8
limits <- aes(ymax = mean.frac + sd.frac, ymin=mean.frac - sd.frac)
m.bar <- ggplot(df.out.lst.meanvar, aes(x = gene.type, y = mean.frac)) + geom_bar(stat = "identity", width = barwidth) + theme_bw() + 
  geom_errorbar(limits, width = barwidth / 2) + theme(aspect.ratio=1.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  
  xlab("") + ylab("Fraction of Genes with Liver-specific DHS peaks")
print(m.bar)
dev.off()

# get ROR counts of each gene



# interesting hits
N.mat.freqs.sub <- subset(N.mat.freqs, (pair == "FOXA2;RORA" | pair == "CUX2;RORA" | pair == "ONECUT1.2;RORA") & motif1 == "atop" & motif2 == "atop" & model == "rhyth")
# N.mat.freqs.sub <- subset(N.mat.freqs, (pair == "FOXA2;RORA" | pair == "CUX2;RORA" | pair == "ONECUT1.2;RORA" | pair == "HNF4A_NR2F1.2;RORA") & motif1 == "atop" & motif2 == "atop" & model == "rhyth")
# N.mat.freqs.sub <- subset(N.mat.freqs,  & motif1 == "atop" & motif2 == "atop" & model == "rhyth")

# get hits by gene names
N.mat.all.sub <- subset(N.mat.all, model == "rhyth" & RORA == "atop" & (ONECUT1.2 == "atop" | CUX2 == "atop" | FOXA2 == "atop"), select = c(gene, peak, model, RORA, ONECUT1.2, CUX2, FOXA2))
N.mat.all.sub <- subset(N.mat.all, model == "rhyth" & RORA == "atop" & (ONECUT1.2 == "atop" | CUX2 == "atop" | FOXA2 == "atop" | AHR_ARNT_ARNT2 == "atop"), select = c(gene, peak, model, RORA, ONECUT1.2, CUX2, FOXA2, AHR_ARNT_ARNT2))
# N.mat.all.sub <- subset(N.mat.all, model == "rhyth" & RORA == "atop" & (ONECUT1.2 == "atop" | CUX2 == "atop" | FOXA2 == "atop" | HNF4A_NR2F1.2 == "atop"), select = c(gene, peak, model, RORA, ONECUT1.2, CUX2, FOXA2, HNF4A_NR2F1.2))
# N.mat.all.sub <- subset(N.mat.all, model == "rhyth" & RORA == "atop", select = c(gene, peak, model))


genes.coop <- as.character(N.mat.all.sub$gene)
genes.clockliver <- as.character(subset(N.mat.all, model == "rhyth")$gene)

# # Top pair ALWAYS ROR?
# motifs.all <- colnames(subset(N.mat.all, select = c(-gene, -peak, -model)))
# top.pairs <- data.frame(motif = motifs.all)
# top.pairs$top.partner <- sapply(top.pairs$motif, function(p){
#   fits.sub <- fits[grepl(p, fits$pair), ] %>% arrange(pval)
#   return(fits.sub$pair[[1]])
# })

# Count top single hits for all motifs
N.mat.long <- melt(N.mat.all, id.vars = c("gene", "peak", "model"), variable.name = "motif")
N.mat.sum <- N.mat.long %>%
  group_by(model, motif) %>%
  filter(value == "atop") %>%
  filter(!duplicated(gene)) %>% 
  summarise(n.genes.with.motif = length(which(value == "atop"))) %>%
  filter(model == "rhyth") %>%
  arrange(desc(n.genes.with.motif))  # ROR is most enriched motif!
N.mat.sum <- OrderDecreasing(N.mat.sum, "motif", "n.genes.with.motif")

# is ROR most significantly enriched by singletons?

N.mat.freqs.sing <- N.mat.long %>%
  group_by(model, motif, value) %>%
  summarise(freq = length(gene))
colnames(N.mat.freqs.sing) <- c("model", "motif", "motif1", "freq")


fits.single <- subset(N.mat.freqs.sing, model != "flat") %>%
  group_by(motif) %>%
  do(RunPoissonModel.sing(.)) %>%
  arrange(pval)

# same result with RORA if we just do cutoff with N.long?? 
sc.cutoff <- 0.5
N.sitecounts <- subset(N.long.filt, gene %in% genes.clockliver) %>%
  group_by(motif2) %>%
  filter(sitecount > sc.cutoff) %>%
  summarise(n.gene = length(gene)) %>%
  arrange(desc(n.gene))


# are they special? amplitudes?
fits.bytiss.livWT <- subset(fits.bytiss, tissue == "Liver_SV129" & gene %in% genes.clockliver)
fits.bytiss.livWT$is.hit <- sapply(as.character(fits.bytiss.livWT$gene), function(g) ifelse(g %in% genes.coop, TRUE, FALSE))
fits.bytiss.livWT <- OrderDecreasing(fits.bytiss.livWT, jfactor = "is.hit", jval = "is.hit")

# fits.bytiss.livWT.sum <- fits.bytiss.livWT %>%
#   group_by(is.hit) %>%

fits.bytiss.livWT$model <- fits.bytiss.livWT$is.hit 
weight.hash <- hash(as.character(subset(fits.long.filt, model == "Liver_SV129")$gene), subset(fits.long.filt, model == "Liver_SV129")$weight)
fits.bytiss.livWT$weight <- sapply(as.character(fits.bytiss.livWT$gene), function(g) weight.hash[[g]])

weight.cutoff <- 0.8
fits.bytiss.livWT$above.weight <- sapply(fits.bytiss.livWT$weight, function(w) ifelse(w >= weight.cutoff, TRUE, FALSE))


elf.target.genes <- as.character(subset(N.mat.all, ELF1.2.4 == "atop" & model == "rhyth")$gene)



# Downstream --------------------------------------------------------------

# Hnf1a is present in both kidney liver, how do we explain
subset(N.mat.freqs.all, pair == "RORA;HNF1A")

foxa2.targs <- GetTargetGenesFromPair(N.mat.all, pair = "FOXA2;RORA")
onecut.targs <- GetTargetGenesFromPair(N.mat.all, pair = "ONECUT1.2;RORA")
cux.targs <- GetTargetGenesFromPair(N.mat.all, pair = "CUX2;RORA")
hnf1a.targs <- GetTargetGenesFromPair(N.mat.all, pair = "RORA;HNF1A")
tead.targs <- GetTargetGenesFromPair(N.mat.all, pair = "RORA;TEAD1")
tfap.targs <- GetTargetGenesFromPair(N.mat.all, pair = "RORA;TFAP2.A.C.")

livtf.targs <- bind_rows(foxa2.targs, onecut.targs, cux.targs)
livtf.targ.genes <- as.character(unique(livtf.targs$gene))

hnf1a.targ.genes <- as.character(unique(hnf1a.targs$gene))

which(hnf1a.targ.genes %in% livtf.targ.genes)  # Hectd2 and Snx10 not in Liv targ genes


subset(N.mat.all, pair == "RORA;HNF1A")

# Is Tead1 and Tfap2 kidney specific rather than liver?


