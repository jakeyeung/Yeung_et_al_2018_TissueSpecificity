# Jake Yeung
# 2015-09-18
# diagnostic_plots.R


# Source ------------------------------------------------------------------

source("scripts/functions/LoadAndHandleData.R")
source("scripts/functions/LoadArray.R")
source('scripts/functions/PlotFunctions.R')

saturation <- function(params, r) params[2] + (params[1] * r) / (params[3] + r)
saturation.inv <- function(params, m) params[3] * (m - params[2]) / (params[2] + params[1] - m)

# Load ---------------------------------------------------------------------

# define dirs
data.dir <- "data"

rna.seq.exprs <- LoadRnaSeq()
array.exprs <- LoadArray(get.norm = TRUE, form = "wide")

# load('data/saturation.fit.select.results.slope.adj.0.7.RERUN.RData', verbose=T)
load("/home/yeung/projects/tissue-specificity/data/RData/sat.fit.select.0.07.RData", verbose = T)

common.genes <- intersect(rownames(array.exprs), rownames(rna.seq.exprs))
common.samples <- intersect(colnames(array.exprs), colnames(rna.seq.exprs))

for (gene in c("Hnf4a")){
  fit.select <- fit.select.list[[gene]]
  fit.used <- fit.select$fit.used  # either saturation or lm
  myfit <- fit.select$myfit
  
  # Get vector of predicted values
  x <- GetFullR(gene, rna.seq.exprs, common.samples)
  y <- unlist(array.exprs[gene, ])
  # x.predict <- seq(min(x)*0.8, max(x)*1.2, length.out=10*length(x))
  y.predict <- seq(min(y) * 0.8, max(y), length.out=10*length(y))
  if (fit.used == "saturation"){
    x.hat <- saturation.inv(coef(myfit), y.predict)
  } else if (fit.used == "lm"){
    x.hat <- linear.inv(coef(myfit), y.predict)
  } else {
    warning("Neither saturation nor lm")
  }
  
  # plot 
  params.str <- paste0(signif(as.vector(coef(myfit)), 2), collapse=",")
  symbols <- GetUnobsObsSymbol(all.samples=colnames(array.exprs), common.samples, unobs=8, obs=1)
  sizes <- GetUnobsObsSymbol(all.samples=colnames(array.exprs), common.samples, unobs=0.25, obs=1)
  plot(x, y, main=paste0("Gene=", gene, " Params=", params.str), pch=symbols, cex=sizes,
       xlab="RNA-Seq DESeq-normalized counts",
       ylab="Microarray normal scale")
  lines(x.hat, y.predict)
  plot(log2(x + 1), log2(y), main=paste0("Gene=", gene, " Params=", params.str), pch=symbols, cex=sizes,
       xlab="RNA-Seq DESeq-normalized counts (log2)",
       ylab="Microarray log2")
  lines(log2(x.hat + 1), log2(y.predict))
  
}