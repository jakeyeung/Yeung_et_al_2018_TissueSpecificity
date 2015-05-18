# PlotGeneAcrossTissues.R

PlotGeneAcrossTissues <- function(dat, jtitle){
  library(ggplot2)
  if (missing(jtitle)){
    jtitle = unique(dat$gene)
  }
  m <- ggplot(dat, aes(x = time, y = exprs,
                       group = experiment, 
                       colour = experiment)) + 
    geom_point() + 
    geom_line() + 
    facet_wrap(~tissue) +
    ggtitle(jtitle) + 
    ylab(label = "log2 expression")
  return(m)
}

PlotTpmAcrossTissues <- function(dat, jtitle, log2.transform=FALSE){
  library(ggplot2)
  if (missing(jtitle)){
    jtitle = unique(dat$gene_name)
  }
  if (log2.transform == FALSE){
    m <- ggplot(dat, aes(x = time, y = tpm)) 
    jylab <- "TPM expression"
  } else {
    m <- ggplot(dat, aes(x = time, y = log2(tpm)))
    jylab <- "log2 TPM expression"
  }
  m <- m + 
    geom_point() + 
    geom_line() + 
    facet_wrap(~tissue) +
    ggtitle(jtitle) + 
    ylab(label = jylab)
  return(m)
}
