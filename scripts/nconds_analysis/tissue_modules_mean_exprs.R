# Jake Yeung
# tissue_modules_mean_exprs.R
# What are mean expression of tissue modules?  
# 2015-11-16

library(dplyr)
library(ggplot2)
library(reshape2)
library(hash)


# Source ------------------------------------------------------------------

source("scripts/functions/PlotGeneAcrossTissues.R")


# Functions ---------------------------------------------------------------

PlotOverlayTimeSeries <- function(dat.long, genes, tissues, jscale = T, jalpha = 0.05, jtitle = ""){
  dat.sub <- subset(dat.long, gene %in% genes & tissue %in% tissues)
  
  # scale and center
  dat.sub <- dat.sub %>%
    group_by(gene, experiment) %>%
    mutate(exprs.scaled = scale(exprs, center = T, scale = jscale))
  
  m <- ggplot(subset(dat.sub), aes(x = time, y = exprs.scaled, group = gene)) + geom_line(alpha = jalpha) + facet_wrap(~experiment) + ggtitle(jtitle)
  m <- m + theme_bw(24) + 
    theme(aspect.ratio=1,
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  if (jscale){
    .ylab <- "Exprs (scaled)"
  } else {
    .ylab <- "Exprs (centered)"
  }
  m <- m + ylab(.ylab) + xlab("Time (CT)")
}

PlotMeanExprsOfModel <- function(dat.mean, genes, jmodel){
  dat.mean.sub <- subset(dat.mean, gene %in% genes)
  m <- ggplot(dat.mean.sub, aes(x = tissue, y = exprs.mean)) + geom_boxplot() + theme_bw(24) + 
    theme(aspect.ratio=1,
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    xlab("") + ylab("Mean expression of gene") +
    ggtitle(paste0("Mean expression level of genes in ", jmodel, " module"))
  return(m)
}

# Load --------------------------------------------------------------------

load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.20.phase_sd_maxdiff_avg.Robj")
load("Robjs/dat.long.fixed_rik_genes.Robj")


# PCA on mean expression --------------------------------------------------

dat.mean <- subset(dat.long, experiment == "rnaseq") %>%
  group_by(gene, tissue, experiment) %>%
  summarise(exprs.mean = mean(exprs))


# Count of genes per model ------------------------------------------------

fits.count <- fits.best %>%
  group_by(model) %>%
  summarise(count = length(gene))
fits.count <- fits.count[order(fits.count$count, decreasing = T), ]


# Take all adrenal-rhythmic genes and plot their rhythms ------------------

adr.genes <- as.character(subset(fits.best, model == "Adr")$gene)

m <- PlotOverlayTimeSeries(dat.long, adr.genes, tissues = "Adr", jalpha = 0.05, jtitle = "Adrenal-specific rhythmic genes")
print(m)

# BFAT --------------------------------------------------------------------

BFAT.genes <- as.character(subset(fits.best, model == "BFAT")$gene)

m <- PlotOverlayTimeSeries(dat.long, BFAT.genes, tissues = "BFAT", jalpha = 0.05, jtitle = "BFAT-specific rhythmic genes")
print(m)


# Mus ---------------------------------------------------------------------

jtiss <- "Mus"
Mus.genes <- as.character(subset(fits.best, model == jtiss)$gene)

m <- PlotOverlayTimeSeries(dat.long, Mus.genes, tissues = jtiss, jscale = T, jalpha = 0.05, jtitle = paste0(jtiss, "-specific rhythmic genes"))
print(m)



# Liver -------------------------------------------------------------------


jtiss <- "Liver"
Mus.genes <- as.character(subset(fits.best, model == jtiss)$gene)

m <- PlotOverlayTimeSeries(dat.long, Mus.genes, tissues = jtiss, jscale = T, jalpha = 0.05, jtitle = paste0(jtiss, "-specific rhythmic genes"))
print(m)


# Is Liver-specific genes liver-specific by mean exprs? -------------------

jtiss <- "Liver"
Liver.genes <- as.character(subset(fits.best, model == jtiss)$gene)

dat.mean.sub <- subset(dat.mean, gene %in% Liver.genes)
dat.mat <- dcast(dat.mean.sub, gene ~ tissue, value.var = "exprs.mean")
rownames(dat.mat) <- dat.mat$gene; dat.mat$gene <- NULL

dat.pca <- prcomp(t(dat.mat), center = T, scale. = T)
screeplot(dat.pca, type = "lines")
barplot(dat.pca$x[, 1])

# boxplot of genes
ggplot(dat.mean.sub, aes(x = tissue, y = exprs.mean)) + geom_boxplot() + theme_bw(24) + 
  theme(aspect.ratio=1,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab("") + ylab("Mean expression of gene") +
  ggtitle(paste0("Mean expression level of genes in ", jtiss, " module"))


# What about Adr BFAT Mus -------------------------------------------------

jtiss <- "Adr"
Adr.genes <- as.character(subset(fits.best, model == jtiss)$gene)
PlotMeanExprsOfModel(dat.mean, Adr.genes, jtiss)

jtiss <- "BFAT"
BFAT.genes <- as.character(subset(fits.best, model == jtiss)$gene)
PlotMeanExprsOfModel(dat.mean, BFAT.genes, jtiss)

jtiss <- "Mus"
Mus.genes <- as.character(subset(fits.best, model == jtiss)$gene)
PlotMeanExprsOfModel(dat.mean, Mus.genes, jtiss)

# Fishers test to say it is not enriched for highly expressed genes -------

plot(density(dat.mean$exprs.mean))
abline(v = 4)

cutoff <- 4  # by eye

# 2 by 2 table: use flat model as background
bg.genes <- as.character(subset(fits.best, model == "")$gene)

# hash to track genes
key <- c(Liver.genes, bg.genes)
val <- c(rep("Liver", length(Liver.genes)), rep("BG", length(bg.genes)))
jhash <- hash(key, val)

dat.mean.sub2 <- subset(dat.mean, gene %in% c(Liver.genes, bg.genes) & tissue != "Liver")

dat.mean.sub2$is.exprs <- sapply(dat.mean.sub2$exprs.mean, function(exprs){
  if (exprs > cutoff){
    return(TRUE)
  } else{
    return(FALSE)
  }
})

dat.mean.sub2$model <- as.factor(sapply(dat.mean.sub2$gene, function(x) jhash[[as.character(x)]]))

liver.exprs.tbl <- table(dat.mean.sub2$is.exprs, dat.mean.sub2$model)

fisher.test(liver.exprs.tbl)


# Is it contamination in Adr? ---------------------------------------------

# lets find genes highly expressed in liver but not expressed in adr

# most varying genes in adrenal gland
dat.var <- subset(dat.long, experiment == "rnaseq" & tissue %in% c("Adr", "Liver")) %>%
  group_by(gene, tissue, experiment) %>%
  summarise(mean.exprs = mean(exprs), var.exprs = var(exprs))

dat.var <- dat.var[order(dat.var$var.exprs, decreasing = TRUE), ]

# plot hits
hits <- dat.var$gene[1:50]
dat.sub <- subset(dat.long, gene %in% hits)



pdf("plots/is_adr_contaminted.pdf")
for (jgene in dat.var$gene){
  print(PlotGeneAcrossTissues(subset(dat.sub, gene == jgene)))
}
dev.off()


# What do rhythmic genes in adrenal gland look like? ----------------------

m <- PlotGeneAcrossTissues(subset(dat.long, gene == "Rgs16"))
m2 <- PlotGeneAcrossTissues(subset(dat.long, gene == "Cps1"))
jlayout <- jlayout <- matrix(c(1, 2, 3, 4), 2, 2, byrow = TRUE)
multiplot(m, m2, layout = jlayout)

dat.mean.adrliv <- subset(dat.mean, tissue %in% c("Liver", "Adr"))

dat.mean.diff.adrliv <- dat.mean.adrliv %>%
  group_by(gene) %>%
  filter(exprs.mean[1] < 0.5 & exprs.mean[2] > 6)

head(dat.mean.diff.adrliv[order(dat.mean.diff.adrliv$exprs.mean, decreasing = T), ])
