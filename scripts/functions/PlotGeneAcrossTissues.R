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
    jgene <- unique(dat$gene_name)
    jtranscript <- unique(dat$transcript_id)
    jtitle = paste(jgene, jtranscript)
  }
  if (log2.transform == FALSE){
    m <- ggplot(dat, aes(x = time, y = tpm, group = transcript_id, colour = transcript_id, shape = transcript_id)) 
    jylab <- "TPM expression"
  } else {
    m <- ggplot(dat, aes(x = time, y = log2(tpm + 0.01), group = transcript_id, colour = transcript_id, shape = transcript_id))
    jylab <- "log2 TPM expression"
  }
  m <- m + theme(legend.position="bottom") +
    geom_point() + 
    geom_line() + 
    facet_wrap(~tissue) +
    ggtitle(jtitle) + 
    ylab(label = jylab)
  return(m)
}

PlotRnaseqAcrossTissues <- function(dat, jtitle){
  p <- ggplot(dat, aes(x = time, y = exprs)) + 
    geom_line() + 
    facet_wrap(~tissue) +
    ggtitle(jtitle) + 
    ylab(label = "log2 mRNA expression") +
    xlab("CT") +
    scale_x_continuous(limits = c(18, 64), breaks = seq(24, 64, 12))  +
    theme(axis.text.x=element_text(angle=90,vjust = 0))
  return(p)
}