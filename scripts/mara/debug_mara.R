# debug_mara.R
# Jake Yeung
# 2015-06-23
# Debug MARA because FOV is too high

library(ggplot2)
# Functions ---------------------------------------------------------------

source("scripts/functions/LoadLong.R")
source("scripts/functions/PlotGeneAcrossTissues.R")

MergedToLong <- function(dat){
  tissues <- sapply(colnames(dat), function(x){
    tisstime <- strsplit(x, "_")[[1]][[1]]
    tiss <- substr(tisstime, start = 1, stop = (nchar(tisstime) - 2))
    return(tiss)
  })
  times <- sapply(colnames(dat), function(x){
    tisstime <- strsplit(x, "_")[[1]][[1]]  
    time <- substr(tisstime, start = nchar(tisstime) - 1, stop = nchar(tisstime))
  })
  experiments <- sapply(colnames(dat), function(x){
    experiment <- strsplit(x, "_")[[1]][[2]]  
  })
  return(data.frame(gene = rep(rownames(dat), ncol(dat)),
                    tissue = rep(tissues, each = nrow(dat)),
                    time = as.numeric(rep(times, each = nrow(dat))),
                    experiment = rep(experiments, each = nrow(dat)),
                    exprs = unlist(as.vector(dat))))
}

# Load --------------------------------------------------------------------

ka.long <- LoadLong()

sitecounts.path <- source("scripts/functions/LoadSitecounts.R")
N.promoter <- LoadSitecounts()

indir <- "/home/yeung/projects/tissue-specificity/results/MARA/rhythmic_genes_by_tissues_kallistoarray.pval0.001.relamp0.1.rerun/"
load(file.path(indir, "Lung", "r.Robj"), verbose = TRUE)  # r, E, N, exp, site
# preds.long <- MergedToLong(preds)
# head(preds.long)

tissues <- dir(indir) 

preds.all.list <- list()
for (tissue in tissues){
  load(file.path(indir, tissue, "r.Robj"), verbose = TRUE)  # r, E, N, exp, site
  preds <- N %*% r$Ahat
  preds.long <- MergedToLong(preds)
  preds.all.list[[tissue]] <- preds.long
}
preds.all <- do.call(what = rbind, args = preds.all.list)


# Plot stuff --------------------------------------------------------------

jgene <- "Arntl"
ggplot(subset(preds.all, gene == jgene), 
       aes(x = time, y = exprs, group = experiment, colour = experiment, fill = experiment)) +
  geom_point() + geom_line() + facet_wrap(~tissue)

PlotGeneAcrossTissues(subset(ka.long, gene == jgene))
