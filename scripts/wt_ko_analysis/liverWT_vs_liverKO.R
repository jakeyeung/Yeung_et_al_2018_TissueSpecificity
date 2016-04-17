# 2016-04-15
# liverWT_vs_liverKO.R
# Do liver-rhythmic genes matter in DD vs LD? Also WT vs KO.

rm(list=ls())

library(ggplot2)
library(dplyr)
library(parallel)

# Load --------------------------------------------------------------------

load("Robjs/dat.long.fixed_rik_genes.Robj", v=T)
load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.15.phase_sd_maxdiff_avg.Robj", v=T)
load("Robjs/dat.fit.Robj", v=T)
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/SvdFunctions.R")
source("scripts/functions/FitRhythmic.R")
source("scripts/functions/NcondsFunctions.R")
source("scripts/functions/Queue.R")
source("scripts/functions/ListFunctions.R")
source("scripts/functions/MemoryManagement.R")
source("scripts/functions/EnvironmentsFunctions.R")

# Read table --------------------------------------------------------------

inf <- "/home/yeung/data/tissue_specificity/cedric_liver_WT_KO/Cry_WT_KO_i_e.txt"
mat.wtko <- read.table(inf, header = TRUE, sep = "\t")

mat.wtko$gene_name <- sapply(as.character(mat.wtko$gene_name), function(g) strsplit(g, "\\|")[[1]][[2]])
genes <- mat.wtko$gene_name
mat.wtko <- mat.wtko[, grepl("_Exon_", colnames(mat.wtko))]
genotype <- sapply(colnames(mat.wtko), function(cname) strsplit(cname, "_")[[1]][[2]])
time <- as.numeric(sapply(colnames(mat.wtko), function(cname) strsplit(cname, "_")[[1]][[5]]))
biorep <- as.numeric(sapply(colnames(mat.wtko), function(cname) strsplit(cname, "_")[[1]][[6]]))
# multiply time by biorep to get the second day 
time.long <- time + (biorep - 1) * 24

dat.wtko <- data.frame(gene = genes,
                       genotype = rep(genotype, each = nrow(mat.wtko)),
                       time = rep(time.long, each = nrow(mat.wtko)),
                       biorep = rep(biorep, each = nrow(mat.wtko)),
                       experiment = "WTKO",
                       tissue = rep(genotype, each = nrow(mat.wtko)), 
                       exprs = unlist(mat.wtko))

fit.wtko <- FitRhythmicDatLong(dat.wtko)

# Define my genes ---------------------------------------------------------

genes.liv <- as.character(subset(fits.best, model == "Liver")$gene)
genes.tw <- as.character(subset(fits.best, n.rhyth >= 8)$gene)

dat.wtko.sub <- subset(dat.wtko, gene %in% genes.liv)

PlotGeneAcrossTissues(subset(dat.wtko, gene == "Dbp"))
PlotGeneAcrossTissues(subset(dat.wtko, gene == "Npas2"))
PlotGeneAcrossTissues(subset(dat.wtko, gene == "Arntl"))
PlotGeneAcrossTissues(subset(dat.wtko, gene == "Cdkn1a"))
PlotGeneAcrossTissues(subset(dat.wtko, gene == "Nr1d1"))
PlotGeneAcrossTissues(subset(dat.wtko, gene == "Saa2"))
PlotGeneAcrossTissues(subset(dat.wtko, gene == "Orm2"))
PlotGeneAcrossTissues(subset(dat.wtko, gene == "Lca5l"))
PlotGeneAcrossTissues(subset(dat.wtko, gene == "Slc2a4"))


# Compare WT and KO in Liver TISSUE WIDE ----------------------------------

dat.test <- TemporalToFrequency2(subset(dat.wtko, gene == "Dbp" & genotype == "WT"), period = 24, n = 12, interval = 4, add.entropy.method = FALSE)
dat.complex.wtko <- TemporalToFrequencyDatLong(dat.wtko, period = 24, n = 12, interval = 4, add.entropy.method = FALSE)

s.tw <- SvdOnComplex(subset(dat.complex.wtko, gene %in% genes.tw), value.var = "exprs.transformed")
eigens.tw <- GetEigens(s.tw, period = 24, comp = 1, label.n = 15, eigenval = TRUE, adj.mag = TRUE, constant.amp = 2)
jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)
multiplot(eigens.tw$u.plot, eigens.tw$v.plot, layout = jlayout)

s.liv <- SvdOnComplex(subset(dat.complex.wtko, gene %in% genes.liv), value.var = "exprs.transformed")
eigens.liv <- GetEigens(s.liv, period = 24, comp = 2, label.n = 15, eigenval = TRUE, adj.mag = TRUE, constant.amp = 2)
jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)
multiplot(eigens.liv$u.plot, eigens.liv$v.plot, layout = jlayout)

ggplot(subset(dat.complex.wtko, gene %in% genes.liv), aes(x = Mod(exprs.transformed) * 2, fill = tissue)) + geom_density(alpha = 0.5)


# Merge with Hogenesch ----------------------------------------------------

# merge with hogenesch then add the new liver WT and KO
dat.wtko.hog <- rbind(subset(dat.wtko, select=c(gene, tissue, time, experiment, exprs)), subset(dat.long, tissue == "Liver"))
dat.wtko.hog$tissue <- factor(as.character(dat.wtko.hog$tissue), levels = unique(dat.wtko.hog$tissue))

print("Chunking data to environment")
dat.env <- DatLongToEnvironment(dat.wtko.hog)

start <- Sys.time()
# Rprof()
fits.all <- mclapply(ls(dat.env), function(gene){
  gene <- "Dbp"
  x <- MakeDesMatRunFitEnv(dat.env, gene, as.character(unique(dat.wtko.hog$tissue)), 
                      n.rhyth.max = 3, w = 2 * pi / 24, 
                      criterion = "BIC", 
                      normalize.weights = TRUE, 
                      cutoff = 1e-5, top.n = NULL, sparse = FALSE)
}, mc.cores = 5)
print(Sys.time() - start)

fits.all.long.wtkohog <- lapply(fits.all, function(x){
  gene <- x$gene
  x$gene <- NULL
  fits.temp.long <- ListToLong(x, gene, top.n = 15)  # TOP N should be MAX models since we dont sort
})
fits.all.long.wtkohog <- do.call(rbind, fits.all.long.wtkohog)

# take min
fits.all.long.wtkohog <- fits.all.long.wtkohog %>%
  group_by(gene) %>%
  filter(weight.raw == min(weight.raw))

print(Sys.time() - start)

outf <- "Robjs/wtko.fits.all.long.Robj"
if (!file.exists(outf)){
  save(fits.all.long.wtkohog, file=outf)
} else {
  print("Skipping saving")
}


# Analyze downstream ------------------------------------------------------

fits.sum <- subset(fits.all.long.wtkohog, model != "") %>%
  group_by(model) %>%
  summarise(n.genes = length(weight))
fits.sum <- fits.sum[order(fits.sum$n.genes, decreasing = TRUE), ]
fits.sum <- OrderDecreasing(fits.sum, "model", jval = "n.genes")
ggplot(fits.sum, aes(x = model, y = n.genes)) + geom_bar(stat = "identity")

# Where are my clock genes? -----------------------------------------------

fits.sum <- subset(fits.all.long.wtkohog, gene %in% genes.tw & model != "") %>%
  group_by(model) %>%
  summarise(n.genes = length(weight))
fits.sum <- fits.sum[order(fits.sum$n.genes, decreasing = TRUE), ]
fits.sum <- OrderDecreasing(fits.sum, "model", jval = "n.genes")
ggplot(fits.sum, aes(x = model, y = n.genes)) + geom_bar(stat = "identity")


# Collapse models ---------------------------------------------------------

# COLLAPSE MODEL
fits.all.long.wtkohog$model <- mapply(FilterModelByAmp, fits.all.long.wtkohog$model, fits.all.long.wtkohog$param.list, MoreArgs = list(amp.cutoff = 0.15))
fits.all.long.wtkohog$n.params <- sapply(fits.all.long.wtkohog$model, function(m) return(length(strsplit(as.character(m), ";")[[1]])))
fits.all.long.wtkohog$n.rhyth <- sapply(fits.all.long.wtkohog$model, GetNrhythFromModel)
fits.all.long.wtkohog$amp.avg <- mapply(GetAvgAmpFromParams, fits.all.long.wtkohog$param.list, fits.all.long.wtkohog$model)
fits.all.long.wtkohog$phase.sd <- mapply(GetSdPhaseFromParams, fits.all.long.wtkohog$param.list, fits.all.long.wtkohog$model)
fits.all.long.wtkohog$phase.maxdiff <- mapply(GetMaxPhaseDiffFromParams, fits.all.long.wtkohog$param.list, fits.all.long.wtkohog$model)
fits.all.long.wtkohog$phase.avg <- mapply(GePhaseFromParams, fits.all.long.wtkohog$param.list, fits.all.long.wtkohog$model)


# Sum again  --------------------------------------------------------------

fits.sum <- subset(fits.all.long.wtkohog, gene %in% genes.tw & model != "") %>%
  group_by(model) %>%
  summarise(n.genes = length(weight))
fits.sum <- fits.sum[order(fits.sum$n.genes, decreasing = TRUE), ]
fits.sum <- OrderDecreasing(fits.sum, "model", jval = "n.genes")
ggplot(fits.sum, aes(x = model, y = n.genes)) + geom_bar(stat = "identity")

fits.sum <- subset(fits.all.long.wtkohog, gene %in% genes.liv & model != "") %>%
  group_by(model) %>%
  summarise(n.genes = length(weight))
fits.sum <- fits.sum[order(fits.sum$n.genes, decreasing = TRUE), ]
fits.sum <- OrderDecreasing(fits.sum, "model", jval = "n.genes")
ggplot(fits.sum, aes(x = model, y = n.genes)) + geom_bar(stat = "identity")

fits.sub.liv <- subset(fits.all.long.wtkohog, gene %in% genes.liv & model != "")

genes.liv.wtliv <- as.character(subset(fits.sub.liv, model %in% c("WT,Liver", "WT;KO;Liver", "WT;Liv"))$gene)
genes.liv.wtliv <- as.character(subset(fits.sub.liv, model %in% c("WT,KO,Liver"))$gene)
# genes.liv.wtliv <- as.character(subset(fits.sub.liv, model == "WT;KO;Liver")$gene)
# genes.liv.wtliv <- as.character(subset(fits.sub.liv, model == "WT;Liv")$gene)

s.liv.wtliv <- SvdOnComplex(subset(dat.complex.wtko, gene %in% genes.liv.wtliv), value.var = "exprs.transformed")
eigens.tw <- GetEigens(s.liv.wtliv, period = 24, comp = 1, label.n = 15, eigenval = TRUE, adj.mag = TRUE, constant.amp = 2)
jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)
multiplot(eigens.tw$u.plot, eigens.tw$v.plot, layout = jlayout)

plot(density(subset(dat.fit, tissue == "Liver" & gene %in% genes.liv.wtliv)$phase))
plot(density(subset(dat.fit, tissue == "Liver" & gene %in% genes.liv)$phase))

