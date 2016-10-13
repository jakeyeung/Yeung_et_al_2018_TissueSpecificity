# PlotGeneAcrossTissues.R

PlotGeneAcrossTissues <- function(dat, jtitle, convert.linear = FALSE, make.pretty = FALSE, jxlab="CT", do.facet.wrap=TRUE, by.linetype=FALSE){
  library(ggplot2)
  if (missing(jtitle)){
    jtitle = unique(dat$gene)
  }
  if (convert.linear){
    dat$exprs <- (2 ^ dat$exprs) - 1 
    jylab <- "mRNA expression (linear scale)"
  } else {
    jylab <- "mRNA expression (log2 scale)"
  }
  n.exper = length(as.character(unique(dat$experiment)))
  
  if (n.exper > 1){
    m <- ggplot(dat, aes(x = time, y = exprs,
                         group = experiment, 
                         colour = experiment)) 
  } else {
    m <- ggplot(dat, aes(x = time, y = exprs)) 
  }
    m <- m + 
      ggtitle(jtitle) + 
      ylab(label = jylab) + 
      xlab(label = jxlab)
    if (do.facet.wrap){
      m <- m + geom_point() + geom_line() + facet_wrap(~tissue)
    } else {
      if (!by.linetype){
        m <- m + geom_point(data = dat, aes(x = time, y = exprs, group = tissue, colour = tissue)) + geom_line(data = dat, aes(x = time, y = exprs, group = tissue, colour = tissue))
      } else {
        m <- m + geom_point(data = dat, aes(x = time, y = exprs, group = tissue, linetype = tissue)) + geom_line(data = dat, aes(x = time, y = exprs, group = tissue, linetype = tissue))
      }
    }
  if (make.pretty){
   m <- m + theme_bw() + theme(panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               aspect.ratio = 1)
  }
  return(m)
}

PlotGeneTissuesWTKO <- function(dat, timelabel="ZT", jtitle="", split.by="geno", ncols = 2){
  # split by geno or tissue
  m <- ggplot(dat, aes(x = time, colour = tissue, linetype = geno, y = exprs)) + 
    geom_point() + geom_line() + xlab(timelabel) + ylab("log2 mRNA accumulation") + 
    theme_bw(24) + ggtitle(jtitle) + 
    theme(aspect.ratio = 1, legend.position = "bottom")
  if (split.by == "geno"){
    m <- m + facet_wrap(~geno, ncol = ncols)
  } else if (split.by == "tissue"){
    m <- m + facet_wrap(~tissue, ncol = ncols)
  } else {
    warning("Split by must be geno or tissue")
  }
  return(m)
}


PlotGeneAcrossTissuesRnaseq <- function(dat, jtitle, convert.linear = FALSE){
  library(ggplot2)
  if (missing(jtitle)){
    jtitle = unique(dat$gene)
  }
  if (convert.linear){
    dat$exprs <- (2 ^ dat$exprs) - 1 
    jylab <- "mRNA expression (linear scale)"
  } else {
    jylab <- "mRNA expression (log2 scale)"
  }
  
  m <- ggplot(dat, aes(x = time, y = exprs)) + 
    geom_line() + 
    facet_wrap(~tissue) +
    ggtitle(jtitle) + 
    ylab(label = jylab) +
    xlab("CT") +
    scale_x_continuous(limits = c(18, 64), breaks = seq(24, 64, 12)) +
    theme(axis.text.x=element_text(angle=90,vjust = 0)) + theme_bw(24)
}

PlotTpmAcrossTissues <- function(dat, jtitle, log2.transform=FALSE, transcript_id = "transcript_id"){
  library(ggplot2)
  
  if (missing(jtitle)){
    jgene <- unique(dat$gene_name)
    jtranscript <- unique(dat$transcript_id)
    jtitle = paste(jgene, jtranscript)
  }
  if (log2.transform == FALSE){
    m <- ggplot(dat, aes_string(x = "time", y = "tpm", group = transcript_id, colour = transcript_id, shape = transcript_id)) 
    jylab <- "TPM expression"
  } else {
    dat$log2tpm <- log2(dat$tpm + 0.01)
    m <- ggplot(dat, aes_string(x = "time", y = "log2tpm", group = transcript_id, colour = transcript_id, shape = transcript_id))
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
    theme_bw(24) + 
    theme(axis.text.x=element_text(angle=90,vjust = 0))
  return(p)
}

Center <- function(x){
  # Running scale on dplyr doesnt work, try it manually
  return(x - mean(x))
}

PlotGeneNormalized <- function(dat, jtitle){
  if (missing(jtitle)){
    jtitle <- unique(dat$gene)
  }
  dat.norm <- dat %>%
    group_by(tissue) %>%
    mutate(exprs.scaled = Center(exprs))
  p <- ggplot(dat.norm, aes(x = time, y = exprs.scaled, group = tissue, colour = tissue, fill = tissue)) + 
    geom_line() +
    geom_point() + 
    ylab("Centered mRNA expression") +
    xlab("CT")
  return(p)
}

PlotEncodeRnaseq <- function(dat, jtitle, sort.by.tissue = TRUE, by.var = "tpm"){
  source("~/projects/tissue-specificity/scripts/functions/SortByTissue.R")
  dat <- SortByTissue(dat, by.var = by.var)
  
  if (missing(jtitle)){
    jtitle <- unique(dat$gene)
  }
  p <- ggplot(dat, aes(x = tissue, y = tpm)) + geom_bar(stat = "identity") + ggtitle(jtitle) + xlab("") + ylab("Gene expression (TPM)") + 
    theme_bw(24) + theme(axis.text.x=element_text(angle=90,vjust = 0, hjust = 1),
                         panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank())
  return(p)
}

CalculatePeriodogramLong <- function(dat, jexperiment = "array", time.interval = NULL, remove.inf = TRUE){
  # interval between sampling points, in hours
  if (jexperiment == "array"){
    interval <- 2  # hrs
  } else if (jexperiment == "rnaseq"){
    interval <- 6  # hrs
  } else {
    interval <- jinterval
  }
  exprs <- dat$exprs
  p <- CalculatePeriodogram(exprs)
  periods <- signif(interval / p$freq, digits = 3)
  dat.var.s <- data.frame(periodogram = p$p.scaled, period = periods)
  dat.var.s$period <- factor(dat.var.s$period, 
                              levels = sort(unique(dat.var.s$period), decreasing = TRUE))
  if (remove.inf){
    dat.var.s <- subset(dat.var.s, period != Inf)
  }
  return(dat.var.s)
}

PlotPeriodogramLong <- function(dat, jexperiment = "array", time.interval = NULL, jtitle = "Plot title"){
  # expect dat to be by gene
  # Plot periodogram of the gene
  source("~/projects/tissue-specificity/scripts/functions/FourierFunctions.R")
  n.exper <- length(as.character(unique(dat$experiment)))
  if (n.exper > 1){
    dat <- subset(dat, experiment == jexperiment)
  }
  
  dat.periodogram <- dat %>%
    group_by(gene, tissue) %>%
    do(CalculatePeriodogramLong(., jexperiment, time.interval))
  ggplot(dat.periodogram, aes(x = period, y = periodogram)) + facet_wrap(~tissue) +  geom_bar(stat = "identity") + 
    theme_bw() + 
    theme(axis.text.x=element_text(angle=90,vjust = 0, hjust = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          aspect.ratio = 1) + 
    xlab("Periodogram") + ylab("Period [h]")
}