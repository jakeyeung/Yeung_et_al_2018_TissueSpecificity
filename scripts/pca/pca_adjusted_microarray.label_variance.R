# from pca_adjusted_microarray.R
# 2015-09-14
# Jake Yeung

library(ggplot2)
library(reshape)
library(wordcloud)  # for showing text without jumbling
library(hash)

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

dat <- LoadNormalizedArray(normalized.array.path = "data/array.adj.0.07.txt", 
                           remove.negs = TRUE, fix.rik.xgene = TRUE)
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

dat <- log(dat, base = 2)

# Calculate PCA and Screeplot ---------------------------------------------

# center data
dat <- t(scale(t(dat), center = TRUE, scale = FALSE))
dat_pca <- prcomp(t(dat), center=FALSE, scale.=FALSE)

# screeplot(dat_pca, type="lines", npcs = min(100, length(dat_pca$sdev)), log="y", main = "")
npcs <- 100
sdev.norm <- sapply(dat_pca$sdev, function(x) x ^ 2 / sum(dat_pca$sdev ^ 2))
plot(x = 1:npcs, 
     sdev.norm[1:npcs], 
     type='o', 
     log = "y", 
     main = paste0("Variance of first ", npcs, " components"),
     xlab = paste0("Components (", length(dat_pca$sdev), " total components)"),
     ylab = "Normalized variance (sum = 1)")


# Consolidate PCA into long -----------------------------------------------

# tissue loadings
jtissues <- GetTissues(rownames(dat_pca$x), get_unique = FALSE)
jtimes <- as.numeric(GetTimes(rownames(dat_pca$x), get_unique = FALSE))
pca.long <- data.frame(tissue = rep(jtissues, times = ncol(dat_pca$x)),
                       time = rep(jtimes, times = ncol(dat_pca$x)),
                       pc = rep(colnames(dat_pca$x), each = nrow(dat_pca$x)),
                       loading = unlist(as.vector(dat_pca$x), use.names = FALSE))
head(pca.long)


# Plot long PCA -----------------------------------------------------------

ggplot(subset(pca.long, pc == "PC1"), aes(x = time, y = loading)) + geom_point() + geom_line() + facet_wrap(~tissue)
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


# Summarize each PCA by its median T.max? ---------------------------------

pca.p.med <- pca.p %>%
  group_by(pc) %>%
  do(SummarisePeriodogram(.))


# Plot screeplot again: but this time label each component by its max period and number of tissues
# label eigen plot
eigenvals <- dat_pca$sdev ^ 2 / sum(dat_pca$sdev ^ 2)
pcs <- c(paste("PC", seq(288), sep = ""))
eigenvals.dic <- hash(pcs, eigenvals)

# adjust factors for plotting
pca.p.med$eigenvals <- sapply(as.character(pca.p.med$pc), function(pc) eigenvals.dic[[pc]])
pca.p.med$pc.num <- sapply(as.character(pca.p.med$pc), function(x) as.numeric(substr(x, 3, nchar(x))))
pca.p.med$pc <- factor(as.character(pca.p.med$pc), levels = pcs)

# discretize T.max.med to either "Inf", 24, 12, or other
pca.p.med$T.max.med <- sapply(pca.p.med$T.max.med, function(x){
  Ts <- c("Inf", "24")
  if (! x %in% Ts){
    x <- "Other"
  }
  return(x)
})

# bin pca.p.med so I can plot many barplots with facet_wrap()
pca.p.med$bini <- BinVector(pca.p.med$pc.num, 10, 100, 288)

# T.max.med should be factor to get colors right
pca.p.med$T.max.med <- factor(pca.p.med$T.max.med, levels = c("Inf", "24", "Other"))

starts <- c(1, 10, 100)
ends <- c(9, 99, 288)

for (i in seq(length(starts))){
  start <- starts[i]
  end <- ends[i]
  break.i <- floor((end - start ) / 4)
  # http://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
  # getting the right colours in different subsets is not trivial!
  m <- ggplot(subset(pca.p.med, pc.num <= end & pc.num >= start), aes(x = pc.num, y = eigenvals, colour = T.max.med, alpha = N / 12)) + 
    geom_bar(stat = "identity", size = 1.1) + scale_colour_discrete(drop=TRUE, limits = levels(pca.p.med$T.max.med)) + scale_x_continuous(breaks=seq(start, end, break.i)) +
    xlab("Component") + ylab("Fraction of total sqr eigenvals")
  print(m)
}

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
