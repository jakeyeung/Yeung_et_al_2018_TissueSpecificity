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

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

PlotPCTissue <- function(dat){
  ggplot(dat, aes(x = time, y = loading)) + geom_point() + geom_line() + facet_wrap(~tissue)
}

SummarisePeriodogram <- function(dat){
  # Summarize periodogram from GetPeriodogramFreq to 
  # get mode T.max, number of tissues with T.max, mean p.pmax for T.max = median(T.max)
  T.max.med <- as.character(Mdoe(dat$T.max))
  dat.sub <- subset(dat, T.max == T.max.med)
  N <- nrow(dat.sub)
  outdat <- data.frame(T.max.med = T.max.med, N = N)
  return(outdat)
}

GetPeriodogramFreq <- function(dat, interval = 2){
  # interval: hour between each sample. Used to get proper period from frequency
  # 2015-09-14
  x <- dat$loading
  p <- CalculatePeriodogram(x, is.matrix = FALSE)
  max.freq <- FindMaxFreqs(freq = p$freq, periodogram = p$p.scaled, n = 1)
  p.pmax <- p$p.scaled[which(p$freq == max.freq)]
  p.24 <- p$p.scaled[which(p$freq == 2/24)]
  p.inf <- p$p.scaled[which(p$freq == 0)]
  p.time <- sum(p$p.scaled[which(p$freq > 0)])
  # return as dataframe so we can use dplyr
  #   outdat <- list(p.24, p.inf, p.pmax, max.freq)
  outdat <- data.frame(p.24 = p.24, p.inf = p.inf, p.time = p.time, p.pmax = p.pmax, T.max = interval / max.freq)
  return(outdat)
}

MeanCenterAcrossGroups <- function(x, n.per.group=24) {
  # Input:
  # x = vector x of inputs across all samples, contains n.per.group
  # elements per group. 
  # Vector is ordered such that all elements
  # in group i are together. 
  # Counting n.per.group would
  # reveal each group. 
  # 
  # Output:
  # x.mean.center <- calculate mean of each group (i.e. mean
  # of every group of n.per.group) and subtract each element
  # in that group by its mean.
  
  N <- length(x)
  x.mean.center <- rep(NA, length(x))
  # BEGIN SANITY CHECK
  # check that N / n.per.group gives no remainder
  if (N %% n.per.group != 0) {
    warning("Uneven number of elements per group.")
  }
  n.groups <- N / n.per.group  # should be an integer
  # END SANITY CHECK
  
  
  for (i in 0:n.groups - 1) {
    # loop from 0 to n.groups - 1 so that our start index
    # loops like 1, 1 + n.per.group, 1 + 2 * n.per.group
    # for example start.index = 1, 25, 49 ...
    start.index <- i * n.per.group + 1
    end.index <- start.index + n.per.group - 1
    # get group e.g. x[1:24]
    group <- x[start.index:end.index]
    # subtract each element in group by mean(group)
    group.mean.center <- group - mean(group)
    # put into x.mean.center
    x.mean.center[start.index:end.index] <- group.mean.center
  }
  
  return(x.mean.center)
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

dat_pca <- prcomp(t(dat), center=TRUE, scale.=TRUE)

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
pcs <- c("PC14", "PC15")
ggplot(subset(pca.long, tissue == "BFAT" & pc %in% pcs), aes(x = time, y = loading)) + geom_point() + geom_line() + facet_wrap(~pc)


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

ggplot(subset(pca.p.med, pc.num < 25 & pc.num >= 9), aes(x = pc.num, y = eigenvals, fill = T.max.med, alpha = N / 12)) + geom_bar(stat = "identity")

# Plot PCA: tissue components ----------------------------------------------------

# Loop to plot scatter plot of PCA i versus PCA i + 1
for (x_comp in 1:2) {
  
  y_comp <- x_comp + 1  # plot x PCA component i on x-axis, PCA component i+1 on y-axis
  
  
  # Show PCA plot of tissues ----------------------------------------------------  
  # 2. 
  # Color by tissue
  colors.by.tissue <- rep(1:12, each=24)
  textplot(dat_pca$x[, x_comp], dat_pca$x[, y_comp], rownames(dat_pca$x), cex=0.7, col=colors.by.tissue, 
           main=paste("Component", x_comp, "vs.", y_comp))
  PlotLoadings(dat_pca$x[, x_comp], title=paste("Vector Loadings for PCA component:", x_comp))
}

# Which vector loadings have oscillating components? ----------------------
# We will fit two linear models:
# 1) linear model with common intercept and common oscillating component
# 2) linear model with tissue-specific intercept 
# and common oscillating component

# Optional: can MeanCenter the vector y to remove cross-tissue differences
# and focus only on oscillations.

N <- nrow(dat_pca$x)  # number of samples.
T <- 24  # 24 hours in a period

# Fit PCA component of interest: loop to try many different PCAs
# user changeable range
for (pca_vector in 10:20) {
  # Create response vector, which is loadings
  
  y <- dat_pca$x[, pca_vector]
  
  # Optional: 
  # y <- MeanCenterAcrossGroups(y)
  
  # BEGIN: plot periodograms to see which frequency has high activity
  freq.and.periodogram <- CalculatePeriodogram(y)  # returns a list
  freq <- freq.and.periodogram$freq
  periodogram <- freq.and.periodogram$p.scaled
  periodogram.unscaled <- freq.and.periodogram$p.unscaled
  
  # Calculate top 5 frequencies
  max.freqs <- FindMaxFreqs(freq, periodogram)
  
  PlotPeriodogram(freq, periodogram, title=paste("Periodogram for PCA component:", pca_vector))
  # add vertical line at max frequency
  max.f <- max.freqs[1]
  # calculate period from frequency
  max.T <- (1 / max.f) * 2  # multiply by 2 because samples are every 2 hours 
  abline(v=max.f, col='blue', lwd=2)
  # add text to show freq and period.
  # x offset of +0.02 so you can see the text
  text(max.f + 0.02, 0, paste0("T=", signif(max.T, digits=3), "hrs"))
  
  PlotLoadings(y, title=paste("Vector Loadings for PCA component:", pca_vector))
  # END: plot periodograms to see which frequency has high activity
  
#   # BEGIN: Linear fit for period of 24 hours
#   # Create sequence of t = [0, 2, 4, ... (N - 1)]
#   t <- seq(0, 2 * N - 1, 2)
#   
#   # set my angular frequency
#   omega <- (2 * pi) / (T)
#   
#   # set my tissue specific factors, matches my response y to a 
#   # specific tissue. 12 tissues, 24 time points
#   n.tissues <- 11  # 11 conditions because WFAT removed
#   n.timepoints <- 24  # 24 time points
#   
#   tissue.factors <- GetTissueSpecificMatrix(n.tissues, n.timepoints)
#   
#   # fit my lm using cos and sin with angular frequency
#   fit <- lm(y ~ cos(omega * t) + sin(omega * t))
#   fit.multi <- lm(y ~ 0 + tissue.factors + cos(omega * t) + sin(omega * t))
#   # END: Linear fit for period of 24 hours
  
#   # Print some statementsy to describe fit
#   cat("*********************************\n")
#   cat(paste0("PCA ", pca_vector, "\n"))
#   cat("\n")
#   
#   cat("Max freqs\n")
#   print(max.freqs)
#   
#   cat("Simple fit:\n")
#   print(anova(fit))
#   cat("\n")
#   
#   cat("Multi fit:\n")
#   print(anova(fit.multi))
#   cat("\n")
#   
#   cat("Comparing simple and multi fit:\n")
#   print(anova(fit, fit.multi))
#   cat("\n")  
}

# Plot combinations of rhythmicity ----------------------------------------

# rhyth.comps <- c(20, 19, 17, 16, 14, 13)  # normal scale
# rhyth.comps <- c(13, 14, 15, 16, 17, 20)  # if log2 
rhyth.comps <- c(14, 20, 16)
library(combinat)
library(PhaseHSV)
tissues <- c("Liver", "Adr", "Kidney", "BFAT", "Aorta", "Mus", "Mus", "Heart")
tissues <- c("Kidney")
combos <- combn(rhyth.comps, 2)
t <- rep(seq(0, 22, 2), 2)  # for one tissue
time.cols <- hsv(PhaseToHsv(2 * pi * t/ 24, 0, 2 *pi), s=1, v=1)
for (tissue in tissues){
  for (i in 1:ncol(combos)){
    combo <- combos[, i]
    x_comp <- combo[1]
    y_comp <- combo[2]
    t <- rep(seq(1, 12), 2)
    plot(dat_pca$x[grepl(tissue, rownames(dat_pca$x)), x_comp], 
         dat_pca$x[grepl(tissue, rownames(dat_pca$x)), y_comp], 
         cex = 0.7, col = "black", 
         main=paste(tissue, "Component", x_comp, "vs", y_comp),
         xlab = paste("Component", x_comp),
         ylab = paste("Component", y_comp))
    text(dat_pca$x[grepl(tissue, rownames(dat_pca$x)), x_comp],
         dat_pca$x[grepl(tissue, rownames(dat_pca$x)), y_comp],
         rownames(dat_pca$x)[which(grepl(tissue, rownames(dat_pca$x)))],
         col = time.cols)
  }  
}


# All tissues -------------------------------------------------------------

# rhyth.comps <- c(20, 19, 17, 16, 14, 13)
# library(combinat)
# tissues <- c("Liver", "Adr", "Kidney", "BFAT", "Aorta", "Mus", "Mus", "Heart")
# tissue <- paste0(tissues, collapse = "|")
# combos <- combn(rhyth.comps, 2)
# t <- rep(rep(seq(0, 22, 2), 2), length(tissues))  # for all tissues
# time.cols <- hsv(PhaseToHsv(2 * pi * t/ 24, 0, 2 *pi), s=1, v=1)
# for (i in 1:ncol(combos)){
#   combo <- combos[, i]
#   x_comp <- combo[1]
#   y_comp <- combo[2]
#   t <- rep(seq(1, 12), 2)
#   plot(dat_pca$x[grepl(tissue, rownames(dat_pca$x)), x_comp], 
#        dat_pca$x[grepl(tissue, rownames(dat_pca$x)), y_comp], 
#        cex = 0.7, col = "black", 
#        main=paste("Component", x_comp, "vs", y_comp),
#        xlab = paste("Component", x_comp),
#        ylab = paste("Component", y_comp))
#   text(dat_pca$x[grepl(tissue, rownames(dat_pca$x)), x_comp],
#        dat_pca$x[grepl(tissue, rownames(dat_pca$x)), y_comp],
#        rownames(dat_pca$x)[which(grepl(tissue, rownames(dat_pca$x)))],
#        col = time.cols)
# }  
