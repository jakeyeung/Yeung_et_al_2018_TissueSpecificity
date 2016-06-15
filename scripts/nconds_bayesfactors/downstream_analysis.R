# 2016-06-10
# Jake Yeung
# downstream_analysis.R

rm(list=ls())

library(dplyr)
library(ggplot2)

source("scripts/functions/ListFunctions.R")
source("scripts/functions/LiverKidneyFunctions.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/NcondsFunctions.R")
source("scripts/functions/SvdFunctions.R")

# Load --------------------------------------------------------------------

indir <- "/home/yeung/projects/tissue-specificity/Robjs/bayes_factors"
inf <- list.files(indir)

fits.long <- expandingList()
for (f in inf){
  method <- strsplit(f, split = "\\.")[[1]][[5]]
  fpath <- file.path("Robjs/bayes_factors", f)
  load(fpath)
  fits.all.long$method <- method
  fits.long$add(fits.all.long)
}
fits.long <- fits.long$as.list()
fits.long <- do.call(rbind, fits.long)

dat <- LoadLivKid()


# Process lowly expressed genes -------------------------------------------

x <- dat$exprs[which(dat$exprs > min(dat$exprs))]
x <- sample(x, 0.01 * length(x))
plot(density(x)); abline(v=1)
source("scripts/functions/MixtureModelFunctions.R")
cutoff <- FindCutoff(x, lambdas = c(0.01, 0.99), mus = c(-1, 4), k = 2, show.fig = TRUE)$maximum

dat.mean <- dat %>%
  group_by(gene) %>%
  summarise(exprs.mean = mean(exprs))

dat.mean.bytiss <- dat %>%
  group_by(gene, tissue) %>%
  summarise(exprs.mean = mean(exprs))

genes.filt <- as.character(subset(dat.mean, exprs.mean < cutoff)$gene)


# Filter lowly expressed genes --------------------------------------------

dat <- subset(dat, ! gene %in% genes.filt)
fits.long <- subset(fits.long, ! gene %in% genes.filt)

# Project to Frequency ----------------------------------------------------

omega <- 2 * pi / 24
dat.freq <- dat %>%
  group_by(gene, tissue) %>%
  do(ProjectToFrequency2(., omega, add.tissue=TRUE))
s <- SvdOnComplex(dat.freq, value.var = "exprs.transformed")
for (i in seq(1)){
  eigens <- GetEigens(s, period = 24, comp = i, label.n = 15, eigenval = TRUE, adj.mag = TRUE, constant.amp = 4, peak.to.trough = TRUE)
  jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)
  multiplot(eigens$u.plot, eigens$v.plot, layout = jlayout)
}

# Summarize ---------------------------------------------------------------

fits.long.filt <- fits.long %>%
  group_by(gene, method) %>%
  filter(weight.raw == min(weight.raw))

PlotGeneAcrossTissues(subset(dat, gene == "Cry1"))


# Add amplitude and phase info --------------------------------------------

fits.long.filt$amp <- sapply(fits.long.filt$param.list, function(p) GetAvgAmpFromParams(params = p, by.model = FALSE))
fits.long.filt$phase <- sapply(fits.long.filt$param.list, function(p) GetPhaseFromParams(params = p, by.model = FALSE))
fits.long.filt$phase.diff <- sapply(fits.long.filt$param.list, function(p) GetMaxPhaseDiffFromParams(params = p, by.model=FALSE))

# Count -------------------------------------------------------------------

fits.sum <- fits.long.filt %>%
  group_by(method, model) %>%
  summarise(count = length(model))

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
ggplot(fits.sum, aes(x = model, fill = method, y = count)) + geom_bar(stat = "identity", position = "dodge", width=0.5) + theme_bw(24) + scale_fill_manual(values=cbPalette)

# Core-clock genes --------------------------------------------------------

jsub <- subset(fits.long.filt, method == "zf")
print(jsub[order(jsub$weight, decreasing = TRUE), ])


# Show kidney-only genes --------------------------------------------------

jsub <- subset(fits.long.filt, method == "zf" & model == "Kidney" & amp > 0.7)
print(jsub[order(jsub$weight, decreasing = TRUE), ])
print(head(data.frame(jsub[order(jsub$amp, decreasing = TRUE), ]), n = 30))
# Itga6 is good
PlotGeneAcrossTissues(subset(dat, gene == "Scarna8"))

# Show Kidney-Liver antiphasic --------------------------------------------

jsub <- subset(fits.long.filt, method == "zf" & model == "Kidney;Liver" & amp > 0.7)
print(jsub[order(jsub$phase.diff, decreasing=TRUE), ])


# Summarize by SVD --------------------------------------------------------

jmeth <- "zf"
genes.tw <- as.character(subset(fits.long.filt, method == jmeth & model %in% c("Kidney;Liver"))$gene)
genes.tw <- as.character(subset(fits.long.filt, method == jmeth & model %in% c("Kidney"))$gene)
genes.tw <- as.character(subset(fits.long.filt, method == jmeth & model %in% c("Liver"))$gene)
jtiss <- "Kidney"
genes.tw <- as.character(subset(fits.long.filt, method == jmeth & model %in% c("Kidney"))$gene)
# filter
genes.tw <- as.character(subset(dat.mean.bytiss, tissue == jtiss & gene %in% genes.tw & exprs.mean > 3)$gene)

s <- SvdOnComplex(subset(dat.freq, gene %in% genes.tw), value.var = "exprs.transformed")
for (i in seq(1)){
  eigens <- GetEigens(s, period = 24, comp = i, label.n = 30, eigenval = TRUE, adj.mag = TRUE, constant.amp = 4, peak.to.trough = TRUE)
  jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)
  multiplot(eigens$u.plot, eigens$v.plot, layout = jlayout)
}

jgene <- "Pcsk9"
jgene <- "Slc34a2"
jgene <- "Ndrg1"
jgene <- "AC159137.1"
jgene <- "Gm11031"
jgene <- "Scarna8"
jgene <- "Gm16550"
jgene <- "Igfbp1"
jgene <- "Mfsd2a"
jgene <- "Ldlr"
jgene <- "Smarca4"
jgene <- "Kank2"
jgene <- "Ell3"
jgene <- "Serinc4"
jgene <- "Gpihbp1"
jgene <- "Mir694"
jgene <- "Igfbp1"
jgene <- "Gm10787"
jgene <- "Slc19a1"
jgene <- "Col18a1"
jgene <- "mt-Tv"
jgene <- "Slc25a25"
jgene <- "Gm16550"
jgene <- "Col13a1"
PlotGeneAcrossTissues(subset(dat, gene == jgene))


# Analyze genes that differ between models --------------------------------



# Take top genes and ask whether they're highly expressed -----------------

top.genes <- head(names(eigens$eigensamp[order(Mod(eigens$eigensamp), decreasing = TRUE)]), n = 100)

pdf("plots/bayes_factors_math_710/kidney_rhythmic_genes.pdf")
for (jgene in top.genes){
  print(PlotGeneAcrossTissues(subset(dat, gene == jgene)) + geom_vline(xintercept=c(8, 20, 32, 44), linetype = "dotted", alpha = 0.25))
}
dev.off()



