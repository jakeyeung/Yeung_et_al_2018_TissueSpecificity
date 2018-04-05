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

library(scales)
mylog_trans <- function(base=exp(1), from=0){
  # from
  # http://stackoverflow.com/questions/22295253/force-bars-to-start-from-a-lower-value-than-0-in-ggplot-geom-bar-in-r
  trans <- function(x) log(x, base)-from
  inv <- function(x) base^(x+from)
  trans_new("mylog", trans, inv, log_breaks(base=base), 
            domain = c(base^from, Inf))
}



# Load data, log transform ------------------------------------------------

load(file = "Robjs/dat.long.fixed_rik_genes.Robj")

log2.transform <- TRUE
outpdf <- "/home/yeung/projects/tissue-specificity/plots/pca/pca_microarray_by_tissue_time.pdf"
dat <- LoadNormalizedArray(normalized.array.path = "data/array.adj.0.07.txt", 
                           remove.negs = TRUE, fix.rik.xgene = TRUE)
# dat <- LoadArray(form = "wide")
# dat <- LoadRnaSeq(handle.duplicates = FALSE)
# load(file = "Robjs/dat.fit.Robj", verbose = T)
# dat.fit.relamp <- GetRelampByGene(dat.fit, by.gene = ref.gene)
# dat.fit.relamp <- dat.fit.relamp %>%
#   group_by(gene) %>%
#   mutate(max.exprs = max(int.rnaseq), min.pval = min(pval), max.relamp = max(relamp), relamp.norm = relamp / sum(relamp), amp.norm = amp / sum(amp))
# genes.exprs <- unique(subset(dat.fit.relamp, max.exprs >= 4)$gene)
# genes.exprs <- FixRikGenes(genes.exprs)
# 
# dat <- as.matrix(dat)


# Remove non-expressed genes ----------------------------------------------

# dat <- dat[which(rownames(dat) %in% genes.exprs), ]

# Optionally log2 transform -----------------------------------------------

if (log2.transform){
  dat <- log(dat + 1, base = 2)
}


# Filter out WFAT ---------------------------------------------------------

dat <- dat[, !grepl("WFAT", colnames(dat))]

dat.lst <- list(dat, dat[, !grepl("BS|Cere|Hypo", colnames(dat))], dat[, grepl("Liver|Kidney", colnames(dat))])

names(dat.lst) <- c("all", "NoNeuralTissues", "LiverKidney")

pca.global.long.lst <- list()

pdf("plots/primetime_plots_full_paper/pca_plots.pdf")

for (jtitle in names(dat.lst)){
  dat <- dat.lst[[jtitle]]
  dat.centered <- dat - apply(dat, 1, mean)
  dat_pca <- prcomp(t(dat.centered), center=FALSE, scale.=FALSE)
  
  
  
  # Consolidate PCA into long -----------------------------------------------
  
  # tissue loadings
  jtissues <- GetTissues(rownames(dat_pca$x), get_unique = FALSE)
  jtimes <- as.numeric(GetTimes(rownames(dat_pca$x), get_unique = FALSE))
  pca.long <- data.frame(tissue = rep(jtissues, times = ncol(dat_pca$x)),
                         time = rep(jtimes, times = ncol(dat_pca$x)),
                         pc = rep(colnames(dat_pca$x), each = nrow(dat_pca$x)),
                         loading = unlist(as.vector(dat_pca$x), use.names = FALSE))
  n.samps <- length(jtissues)
  
  
  # Analyze whether each PC is rhythmic or not ------------------------------
  
  pca.p <- pca.long %>%
    group_by(pc, tissue) %>%
    do(GetPeriodogramFreq(.))
  
  
  
  # Summarize by number of tissues with 12, 24, Inf, or other ---------------
  
  pca.p.sum <- pca.p %>%
    group_by(pc) %>%
    do(SummarisePeriodogram2(., weighted = TRUE))
  
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
  
  
  # Global tissue differences and temporal differences ----------------------
  
  pca.global <- with(pca.p.sum, data.frame(T.inf.global = sum(T.inf * eigenvals), 
                                           T.12.global = sum(T.12 * eigenvals), 
                                           T.24.global = sum(T.24 *eigenvals), 
                                           T.other = sum(T.other * eigenvals)))
  
  
  # Plot global contributions -----------------------------------------------
  
  pca.global.long <- melt(pca.global, variable.name = "Fourier", value.name = "Fraction")
  pca.global.long <- subset(pca.global.long, Fourier != "T.other")
  pca.global.long$Fourier <- factor(as.character(pca.global.long$Fourier), levels = c("T.inf.global", "T.24.global", "T.12.global"))
  
  pca.global.long$ds <- jtitle
  
  pca.global.long.lst[[jtitle]] <- pca.global.long
  
  print(ggplot(pca.global.long, aes(x = Fourier, y = Fraction)) + 
          geom_bar(stat = "identity") + 
          theme_bw(24) + 
          theme(aspect.ratio = 1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
          ggtitle(jtitle) + 
          ylim(c(0, 1)))
  print(ggplot(subset(pca.global.long, Fourier != "T.inf.global"), aes(x = Fourier, y = Fraction)) + 
          geom_bar(stat = "identity") + 
          theme_bw(24) + 
          theme(aspect.ratio = 1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
          ggtitle(jtitle))
}

print(pca.global.long.lst)

pca.global.long <- do.call(rbind, pca.global.long.lst)

pca.global.long$ds <- factor(pca.global.long$ds, levels = c("all", "NoNeuralTissues", "LiverKidney"))

cols <- c("grey10", "grey50", "grey60")

print(ggplot(pca.global.long, aes(x = Fourier, fill = ds, y = Fraction)) + 
        geom_bar(stat = "identity", position="dodge") +
        theme_bw(24) + 
        theme(aspect.ratio = 1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + 
        scale_y_continuous(trans = mylog_trans(base=10, from=-3), breaks = c(0.001, 0.01, 0.1, 1)) + 
        scale_fill_manual(values = cols))
print(ggplot(subset(pca.global.long, Fourier != "T.12.global"), aes(x = Fourier, fill = ds, y = Fraction)) + 
        geom_bar(stat = "identity", position="dodge") +
        theme_bw(24) + 
        theme(aspect.ratio = 1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + 
        scale_y_continuous(trans = mylog_trans(base=10, from=-2), breaks = c(0.01, 0.1, 1)) + 
        scale_x_discrete(labels=c("T.inf.global" = "Intertissue", "T.24.global" = "Intratissue\n(24h Period)")) +  
        xlab("") +  ylab("Fraction of Total Variance") + 
        scale_fill_manual(values = cols))
print(ggplot(subset(pca.global.long, Fourier != "T.12.global"), aes(x = Fourier, fill = ds, y = Fraction)) + 
        geom_bar(stat = "identity", position="dodge") +
        theme_bw(24) + 
        theme(aspect.ratio = 1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
        scale_y_continuous(trans = mylog_trans(base=10, from=-2), breaks = c(0.01, 0.1, 1)) + 
        scale_x_discrete(labels=c("T.inf.global" = "Intertissue", "T.24.global" = "Intratissue\n(24h Period)")) +  
        xlab("") +  ylab("Fraction of Total Variance") + 
        scale_fill_manual(values = cols))

print(ggplot(pca.global.long, aes(x = Fourier, fill = ds, y = Fraction)) + 
        geom_bar(stat = "identity", position="dodge") +
        theme_bw(24) + 
        theme(aspect.ratio = 1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + 
        ylim(c(0, 1)) + 
        scale_fill_manual(values = cols))
print(ggplot(pca.global.long, aes(x = Fourier, fill = ds, y = Fraction)) + 
        geom_bar(stat = "identity", position="dodge") +
        theme_bw(24) + 
        theme(aspect.ratio = 0.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + 
        ylim(c(0, 1)) + 
        scale_fill_manual(values = cols))
print(ggplot(pca.global.long, aes(x = Fourier, fill = ds, y = Fraction)) + 
        geom_bar(stat = "identity", position="dodge") +
        theme_bw(24) + 
        theme(aspect.ratio = 0.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
        ylim(c(0, 1)) + 
        scale_fill_manual(values = cols))
dev.off()




