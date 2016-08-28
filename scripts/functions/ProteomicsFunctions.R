LoadProteomicsData <- function(inf = "/home/shared/nuclear_proteomics/nuclear_proteins_L_H_log2_all_WT_KO_24h_12h_statistics.OneGenePerLine.txt", 
                               as.long = TRUE){
  require(reshape2)
  require(dplyr)
  prot <- read.table(inf, header = TRUE, sep = "\t")
  if (!as.long){
    return(prot)
  }
  # Make long 
  wt.sampnames <- paste("ZT", sprintf("%02d", seq(0, 45, 3)), ".WT", sep = "")
  ko.sampnames <- paste("ZT", sprintf("%02d", seq(0, 18, 6)), ".Bmal.WT", sep = "")
  fit.sampnames <- c("mean", "amp", "relamp", "phase", "pval", "qv", "amp.12h", "relamp.12h", "phase.12h", "pval.12h", "qv.12h")
  
  prot.long.wt <- melt(prot, id.vars = "Gene.names", measure.vars = wt.sampnames, variable.name = "samp", value.name = "rel.abund")
  prot.long.bmalko <- melt(prot, id.vars = "Gene.names", measure.vars = ko.sampnames, variable.name = "samp", value.name = "rel.abund")
  fit.prot.wt <- subset(prot, select = c("Gene.names", fit.sampnames))
  # fit.prot.wt <- melt(prot, id.vars = "Gene.names", measure.vars = fit.sampnames)
  
  prot.long.wt$time <- GetTimeFromSamp(as.character(prot.long.wt$samp))
  prot.long.wt$geno <- GetGenoFromSamp(as.character(prot.long.wt$samp))
  
  prot.long.bmalko$time <- GetTimeFromSamp(as.character(prot.long.bmalko$samp))
  prot.long.bmalko$geno <- GetGenoFromSamp(as.character(prot.long.bmalko$samp))
  
  # merge
  prot.long <- rbind(prot.long.wt, prot.long.bmalko)
  
  # change Gene.names to gene
  colnames(fit.prot.wt)[which(colnames(fit.prot.wt) == "Gene.names")] <- "gene"
  colnames(prot.long)[which(colnames(prot.long) == "Gene.names")] <- "gene"
  return(prot.long)
}

PlotProteomics <- function(prot.long, jtitle = ""){
  # plot proteomics
  g <- ggplot(prot.long, aes(x = time, y = rel.abund, colour = geno, group = geno)) + geom_point() + geom_line() + theme_bw() 
  g <- g + ggtitle(jtitle)
  return(g)
}

GetTimeFromSamp <- Vectorize(function(samp){
  time <- strsplit(samp, "\\.")[[1]][[1]]
  # remove ZT
  time <- gsub("ZT", "", time)
  return(as.numeric(time))
}, vectorize.args = "samp")

GetGenoFromSamp <- Vectorize(function(samp){
  geno <- strsplit(samp, "\\.")[[1]][[2]]
  return(geno)
}, vectorize.args = "samp")

PlotmRNAActivityProtein <- function(dat.long, act.long, prot.long, gene.dat, gene.act, gene.prot){
  dat.sub <- subset(dat.long, gene == gene.dat & tissue == "Liver")
  act.sub <- subset(act.long, gene == gene.act & tissue %in% c("Liver_SV129", "Liver_BmalKO"))
  prot.sub <- subset(prot.long, gene == gene.prot)
  
  # scale data, merge together, then plot in one figure
  dat.sub$signal <- ScaleSignal(dat.sub, "exprs")
  act.sub$signal <- ScaleSignal(act.sub, "exprs")
  prot.sub$signal <- ScaleSignal(prot.sub, "rel.abund")
  
  # make two factors: WT and BmalKO
  dat.sub$geno.std <- as.factor(gsub("SV129", "WT", as.character(dat.sub$geno)))
  act.sub$geno.std <- as.factor(gsub("SV129", "WT", as.character(act.sub$geno)))
  prot.sub$geno.std <- as.factor(gsub("Bmal", "BmalKO", as.character(prot.sub$geno)))
  
  dat.sub$type <- "mRNA_Accum"
  act.sub$type <- "TF_Activity"
  prot.sub$type <- "Nuclear_Prot_Accum"
  
  dat.sub2 <- subset(dat.sub, select = c(geno.std, signal, time, type))
  act.sub2 <- subset(act.sub, select = c(geno.std, signal, time, type))
  prot.sub2 <- subset(prot.sub, select = c(geno.std, signal, time, type))
  merged.dat <- rbind(dat.sub2,
                      act.sub2,
                      prot.sub2)
  merged.dat$geno.std <- factor(as.character(merged.dat$geno.std), levels = c("WT", "BmalKO"))
  merged.dat$type <- as.factor(merged.dat$type)
  m <- ggplot(merged.dat, aes(x = time, y = signal, linetype = type, shape = type)) + 
    geom_line() +
    geom_point(size = 3) + 
    xlab("Time (ZT)") + ylab("Scaled Signal") + 
    facet_wrap(~geno.std) + theme_bw(24) + 
    theme(legend.position = "bottom", aspect.ratio = 1)
  return(m)
}

ScaleSignal <- function(dat, cname){
  return(scale(dat[[cname]], center = TRUE, scale = TRUE))
}
