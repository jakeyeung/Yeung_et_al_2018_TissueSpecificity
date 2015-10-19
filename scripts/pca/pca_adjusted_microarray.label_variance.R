# from pca_adjusted_microarray.R
# 2015-09-14
# Jake Yeung

rm(list=ls())

library(ggplot2)
library(hash)
library(dplyr)

ref.gene <- "Nr1d1"
# Define my functions -----------------------------------------------------

# First source my functions I wrote
funcs.dir <- file.path('scripts', 'functions')
source(file.path(funcs.dir, 'SampleNameHandler.R'))  # for shortening sample names
source(file.path(funcs.dir, 'PcaPlotFunctions.R'))  # for visualizing PCA, periodoigrams
source(file.path(funcs.dir, 'FourierFunctions.R'))  # for periodoigrams
source(file.path(funcs.dir, 'GetTissueSpecificMatrix.R'))  # as name says
source(file.path(funcs.dir, "GrepRikGenes.R"))
source(file.path(funcs.dir, "GetTissueTimes.R"))
source(file.path(funcs.dir, "PCAFunctions.R"))
source(file.path(funcs.dir, "LoadAndHandleData.R"))
source(file.path(funcs.dir, "FitRhythmic.R"))
source(file.path(funcs.dir, "PlotGeneAcrossTissues.R"))
source(file.path(funcs.dir, "LoadArray.R"))

BinVector <- function(x, ...){
  # bin vector into discrete bins i
  ends <- c(...)
  starts <- c(1, ends[1:length(ends) - 1])
  bins <- sapply(x, function(y){
    for (i in seq(length(starts))){
      if (y >= starts[i] & y <= ends[i]){
        return(i)
      } else {
        next
      }
    }
  })
  return(unlist(bins))
}


# Load data, log transform ------------------------------------------------

load(file = "Robjs/dat.long.fixed_rik_genes.Robj")

log2.transform <- TRUE
outpdf <- "/home/yeung/projects/tissue-specificity/plots/pca/pca_microarray_by_tissue_time.pdf"
dat <- LoadNormalizedArray(normalized.array.path = "data/array.adj.0.07.txt", 
                           remove.negs = TRUE, fix.rik.xgene = TRUE)
# dat <- LoadArray(form = "wide")
# dat <- LoadRnaSeq(handle.duplicates = FALSE)
load(file = "Robjs/dat.fit.Robj", verbose = T)
dat.fit.relamp <- GetRelampByGene(dat.fit, by.gene = ref.gene)
dat.fit.relamp <- dat.fit.relamp %>%
  group_by(gene) %>%
  mutate(max.exprs = max(int.rnaseq), min.pval = min(pval), max.relamp = max(relamp), relamp.norm = relamp / sum(relamp), amp.norm = amp / sum(amp))
genes.exprs <- unique(subset(dat.fit.relamp, max.exprs >= 4)$gene)
genes.exprs <- FixRikGenes(genes.exprs)

dat <- as.matrix(dat)


# Remove non-expressed genes ----------------------------------------------

dat <- dat[which(rownames(dat) %in% genes.exprs), ]

# Optionally log2 transform -----------------------------------------------

if (log2.transform){
  dat <- log(dat + 1, base = 2)
}


# Filter out WFAT ---------------------------------------------------------

dat <- dat[, !grepl("WFAT", colnames(dat))]

# Calculate PCA and Screeplot ---------------------------------------------

# center data
# dat <- t(scale(t(dat), center = TRUE, scale = FALSE))
# center by ROW
dat.centered <- dat - apply(dat, 1, mean)
dat_pca <- prcomp(t(dat.centered), center=FALSE, scale.=FALSE)


# screeplot(dat_pca, type="lines", npcs = min(100, length(dat_pca$sdev)), log="y", main = "")
npcs <- 23
sdev.norm <- sapply(dat_pca$sdev, function(x) x ^ 2 / sum(dat_pca$sdev ^ 2))
plot(x = 1:npcs, 
     sdev.norm[1:npcs], 
     type='o', 
     # log = "y", 
     main = paste0("Variance of first ", npcs, " components"),
     xlab = paste0("Components (", length(dat_pca$sdev), " total components)"),
     ylab = "Normalized variance (sum = 1)")


# PCA1 and PCA2 -----------------------------------------------------------

library(wordcloud)
# plot(dat_pca$x[, "PC1"], dat_pca$x[, "PC2"])
textplot(x = dat_pca$x[, "PC1"], y = dat_pca$x[, "PC2"], words = rownames(dat_pca$x))
# textplot(x = dat_pca$rotation[, "PC17"], y = dat_pca$rotation[, "PC18"], words = rownames(dat_pca$rotation))


# PCA 17 has outlier in Adrenal observe the gene that is responsib --------


jpc <- "PC13"
sorted.eigensamp <- dat_pca$rotation[, jpc][order(abs(dat_pca$rotation[, jpc]), decreasing = TRUE)]
sorted.eigengene <- dat_pca$x[, jpc][order(abs(dat_pca$x[, jpc]), decreasing = TRUE)]
head(sorted.eigengene)
plot(dat_pca$x[, jpc])

maxi <- 30
dat.sub <- subset(dat.long, gene %in% names(sorted.eigensamp[1:maxi]))
pdf(paste0("plots/pca/pc_", jpc, "genes.pdf"))

for (i in seq(30)){
  print(PlotGeneAcrossTissues(subset(dat.sub, gene == names(sorted.eigensamp)[[i]])))
}
dev.off()




# Consolidate PCA into long -----------------------------------------------

# tissue loadings
jtissues <- GetTissues(rownames(dat_pca$x), get_unique = FALSE)
jtimes <- as.numeric(GetTimes(rownames(dat_pca$x), get_unique = FALSE))
pca.long <- data.frame(tissue = rep(jtissues, times = ncol(dat_pca$x)),
                       time = rep(jtimes, times = ncol(dat_pca$x)),
                       pc = rep(colnames(dat_pca$x), each = nrow(dat_pca$x)),
                       loading = unlist(as.vector(dat_pca$x), use.names = FALSE))
head(pca.long)
n.samps <- length(jtissues)

# Plot long PCA -----------------------------------------------------------

# ggplot(subset(pca.long, pc == "PC1"), aes(x = time, y = loading)) + geom_point() + geom_line() + facet_wrap(~tissue)
# ggplot(subset(pca.long, pc == "PC2"), aes(x = time, y = loading)) + geom_point() + geom_line() + facet_wrap(~tissue)
# ggplot(subset(pca.long, pc == "PC3"), aes(x = time, y = loading)) + geom_point() + geom_line() + facet_wrap(~tissue)
# ggplot(subset(pca.long, pc == "PC4"), aes(x = time, y = loading)) + geom_point() + geom_line() + facet_wrap(~tissue)
ggplot(subset(pca.long, pc == "PC17"), aes(x = time, y = loading)) + geom_point() + geom_line() + facet_wrap(~tissue)

pcs <- seq(50); pcs <- paste("PC", pcs, sep = "")

pdf(outpdf)
for (jpc in pcs){
  print(ggplot(subset(pca.long, pc == jpc), aes(x = time, y = loading)) + geom_point() + geom_line() + facet_wrap(~tissue) + ggtitle(jpc))
}
dev.off()

ggplot(subset(pca.long, pc == "PC14" & tissue != "WFAT"), aes(x = time, y = loading)) + geom_point() + geom_line() + facet_wrap(~tissue)
ggplot(subset(pca.long, pc == "PC15"), aes(x = time, y = loading)) + geom_point() + geom_line() + facet_wrap(~tissue)
pcs <- c("PC13", "PC14", "PC15", "PC16", "PC17", "PC18", "PC19", "PC20", "PC21")
ggplot(subset(pca.long, tissue == "Liver" & pc %in% pcs), aes(x = time, y = loading)) + geom_point() + geom_line() + facet_wrap(~pc)
pcs.pair <- c("PC20", "PC21")

x <- subset(pca.long, tissue == "Liver" & pc == pcs.pair[1])$loading
y <- subset(pca.long, tissue == "Liver" & pc == pcs.pair[2])$loading
time <- subset(pca.long, tissue == "Liver" & pc == pcs.pair[2])$time
dat.plot <- data.frame(x = x, y = y, time = time %% 24)
ggplot(dat.plot, aes(x = x, y = y, label = time)) + geom_text()


# Analyze whether each PC is rhythmic or not ------------------------------

pca.p <- pca.long %>%
  group_by(pc, tissue) %>%
  do(GetPeriodogramFreq(.))



# Summarize by number of tissues with 12, 24, Inf, or other ---------------

pca.p.sum <- pca.p %>%
  group_by(pc) %>%
  do(SummarisePeriodogram2(., weighted = TRUE))
# subset(pca.p.sum, pc == "PC13")
# sum(subset(pca.p.sum, pc == "PC13")[1, 2:5])

pca.p.sum$pc.num <- sapply(pca.p.sum$pc, function(p) as.numeric(strsplit(as.character(p), split = "PC")[[1]][[2]]))
pca.p.sum <- pca.p.sum[order(pca.p.sum$pc.num), ]
head(data.frame(pca.p.sum), n = 25)

# Plot screeplot, label each component by number of tissues in eac --------

eigenvals <- dat_pca$sdev ^ 2 / sum(dat_pca$sdev ^ 2)
pcs <- c(paste("PC", seq(n.samps), sep = ""))
eigenvals.dic <- hash(pcs, eigenvals)

# adjust factors for plotting
pca.p.sum$eigenvals <- sapply(as.character(pca.p.sum$pc), function(pc) eigenvals.dic[[pc]])
pca.p.sum$pc.num <- sapply(as.character(pca.p.sum$pc), function(x) as.numeric(substr(x, 3, nchar(x))))
pca.p.sum$pc <- factor(as.character(pca.p.sum$pc), levels = pcs)




# Make long and plot -----------------------------------------------------

pca.p.sum.long <- melt(pca.p.sum, id.var = c("pc", "eigenvals", "pc.num"), value.name = "fracFourier", variable.name = "Component")

pca.p.sum.long <- pca.p.sum.long %>%
  mutate(eigenvals.frac = eigenvals * fracFourier)

starts <- c(1, 13, 100)
ends <- c(12, 99, 288)

for (i in seq(length(starts))){
  print(PlotComponents2(pca.p.sum.long, starts[[i]], ends[[i]]))
}


# Parseval's Theorem ------------------------------------------------------

pca.p.parseval <- pca.p %>%
  group_by(pc) %>%
  summarise(eigen.est = sum(p.inf + p.time) / length(p.inf))
pca.p.parseval$pc.num <- sapply(pca.p.parseval$pc, function(p) as.numeric(strsplit(as.character(p), "PC")[[1]][[2]]))
pca.p.parseval <- pca.p.parseval[order(pca.p.parseval$pc.num), ]

eigenvals <- dat_pca$sdev ^ 2

diff.norm <- abs(eigenvals - pca.p.parseval$eigen.est) ^ 2 / eigenvals

plot(density(diff.norm))

plot(pca.p.parseval$eigen.est, eigenvals)

# # Summarize each PCA by its median T.max? ---------------------------------
# 
# pca.p.med <- pca.p %>%
#   group_by(pc) %>%
#   do(SummarisePeriodogram(.))
# 
# 
# # Plot screeplot again: but this time label each component by its max period and number of tissues
# # label eigen plot
# 
# eigenvals <- dat_pca$sdev ^ 2 / sum(dat_pca$sdev ^ 2)
# pcs <- c(paste("PC", seq(n.samps), sep = ""))
# eigenvals.dic <- hash(pcs, eigenvals)
# 
# # adjust factors for plotting
# pca.p.med$eigenvals <- sapply(as.character(pca.p.med$pc), function(pc) eigenvals.dic[[pc]])
# pca.p.med$pc.num <- sapply(as.character(pca.p.med$pc), function(x) as.numeric(substr(x, 3, nchar(x))))
# pca.p.med$pc <- factor(as.character(pca.p.med$pc), levels = pcs)
# 
# # discretize T.max.med to either "Inf", 24, 12, or other
# pca.p.med$T.max.med <- sapply(pca.p.med$T.max.med, function(x){
#   Ts <- c("Inf", "24")
#   if (! x %in% Ts){
#     x <- "Other"
#   }
#   return(x)
# })
# 
# # bin pca.p.med so I can plot many barplots with facet_wrap()
# pca.p.med$bini <- BinVector(pca.p.med$pc.num, 10, 100, 288)
# 
# # T.max.med should be factor to get colors right
# pca.p.med$T.max.med <- factor(pca.p.med$T.max.med, levels = c("Inf", "24", "Other"))
# 
# starts <- c(1, 13, 100)
# ends <- c(12, 99, 288)
# 
# for (i in seq(length(starts))){
#   PlotComponents(pca.p.med, starts[i], ends[i])
# }

# # replacte gray with some dark blue
# # http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/#a-colorblind-friendly-palette
# cbPalette <- c("#CC79A7", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# for (i in seq(length(starts))){
#   start <- starts[i]
#   end <- ends[i]
#   break.i <- floor((end - start ) / 4)
#   # http://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
#   # getting the right colours in different subsets is not trivial!
#   m <- ggplot(subset(pca.p.med, pc.num <= end & pc.num >= start), aes(x = pc.num, y = eigenvals, colour = T.max.med, alpha = N / 12)) + 
#     geom_bar(stat = "identity", size = 1.1) + scale_colour_discrete(drop=TRUE, limits = levels(pca.p.med$T.max.med)) + scale_x_continuous(breaks=seq(start, end, break.i)) +
#     xlab("Component") + ylab("Fraction of total sqr eigenvals") +
#     scale_colour_manual(values=cbPalette)
#   print(m)
# }

# # ggplot(pca.p.med, aes(x = pc, y = eigenvals, fill = T.max.med, alpha = N / 12)) + geom_bar(stat = "identity") + facet_wrap(~bini)
# 
# 
# ggplot(subset(pca.p.med, pc.num < 25 & pc.num >= 9), aes(x = pc.num, y = eigenvals, colour = T.max.med, alpha = N / 12)) + geom_bar(stat = "identity")
# ggplot(subset(pca.p.med, pc.num < 288 & pc.num >= 100), aes(x = pc.num, y = eigenvals, colour = T.max.med, alpha = N / 12)) + geom_bar(stat = "identity")
# ggplot(subset(pca.p.med, pc.num < 100 & pc.num >= 10), aes(x = pc.num, y = eigenvals, fill = T.max.med, alpha = N / 12)) + geom_bar(stat = "identity")
# ggplot(subset(pca.p.med, pc.num < 288 & pc.num >= 25), aes(x = pc.num, y = eigenvals, fill = T.max.med, alpha = N / 12)) + geom_bar(stat = "identity")
# ggplot(subset(pca.p.med, pc.num < 100 & pc.num >= 10), aes(x = pc.num, y = eigenvals, fill = T.max.med, alpha = N / 12)) + geom_bar(stat = "identity")
# ggplot(subset(pca.p.med, pc.num < 10 & pc.num >= 1), aes(x = pc.num, y = eigenvals, fill = T.max.med, alpha = N / 12)) + geom_bar(stat = "identity")
# ggplot(subset(pca.p.med, pc.num < 25), aes(x = pc.num, y = eigenvals, fill = T.max.med, alpha = N / 12)) + geom_bar(stat = "identity")
