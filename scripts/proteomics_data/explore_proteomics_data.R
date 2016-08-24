# 2016-08-16
# Jake Yeung
# explore_proteomics_data.R

rm(list=ls())

library(dplyr)
library(ggplot2)
library(reshape2)


# Functions ---------------------------------------------------------------

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

# Load --------------------------------------------------------------------

inf <- "/home/shared/nuclear_proteomics/nuclear_proteins_L_H_log2_all_WT_KO_24h_12h_statistics.OneGenePerLine.txt"

prot <- read.table(inf, header = TRUE, sep = "\t")


# Make long ---------------------------------------------------------------

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

# Explore rhythms in proteomics data --------------------------------------

jgene <- "Stat2"
jgene <- "Srfbp1"
jgene <- "Nfil3"
jgene <- "Nr3c1"
jgene <- "Nr3c2"
jgene <- "Atf2"
jgene <- "Tfdp1"
jgene <- "Hmga2"
PlotProteomics(subset(prot.long, gene == jgene)) + ggtitle(jgene)

jgenes <- c("Stat2", "Mafb", "Jun", "E2f5", "Creb3", "Elf2", "Tgif1", "Egr1", "Pou2f1", "Hes6", "Tfdp1", "Crebbp")

jgene <- "Onecut1"
jgene <- "Onecut2"
jgene <- "Nfya"
jgene <- "Nfyc"
jgene <- "Hic1"
# print(PlotProteomics(jsub, jtitle = jgene))

for (jgene in jgenes){
  print(jgene)
  jsub <- subset(prot.long, gene == jgene)
  if (nrow(jsub) > 0 & !all(is.na(jsub$rel.abund))){
    print(PlotProteomics(jsub, jtitle = jgene))
  } else {
    print(paste("Skipping", jgene))
  }
}
