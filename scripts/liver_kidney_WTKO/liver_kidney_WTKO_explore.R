# 2016-06-22
# Jake Yeung
# liver_kidney_WTKO_explore.R

rm(list=ls())

library(dplyr)
library(ggplot2)

source("scripts/functions/BiomartFunctions.R")
source("scripts/functions/PlotGeneAcrossTissues.R")

eps <- 1  # for log2 transform

# Load --------------------------------------------------------------------

inf="/home/shared/atgerWTKO_kidneyWTKO/LiverKidney_SV129BmalKO_RF_TotalRNAPolyARNA_kallisto_abundances.txt"
dat <- read.table(inf, header = TRUE, row.names = 1)

genes <- Transcript2Gene(rownames(dat), return.original = FALSE)  # slow
sampnames <- colnames(dat)
dat$gene <- genes

# Convert to long ---------------------------------------------------------

tissues <- sapply(sampnames, function(s) strsplit(s, "_")[[1]][[1]], USE.NAMES = FALSE)
experiments <- sapply(sampnames, function(s) strsplit(s, "_")[[1]][[2]], USE.NAMES = FALSE)
times <- sapply(sampnames, function(s) strsplit(s, "_")[[1]][[3]], USE.NAMES = FALSE)
genos <- sapply(sampnames, function(s) strsplit(s, "_")[[1]][[4]], USE.NAMES = FALSE)
feedings <- sapply(sampnames, function(s){
  jrep <- strsplit(s, "_")[[1]][[5]]
  # remove .A
  jrep <- strsplit(jrep, "\\.")[[1]][[1]]
  return(jrep)
}, USE.NAMES = FALSE)
replicates <- sapply(sampnames, function(s){
  r <- tryCatch({
    rep.letter <- strsplit(s, "\\.")[[1]][[2]]
    rep.numb <- as.numeric(chartr("ABCDE", "12345", rep.letter))
    return(rep.numb)
  }, error = function(e) {
    return(NA)
  })
  return(r)
}, USE.NAMES = FALSE)

# Integrate replicates into ZT time
times.new <- mapply(function(time, jrep){
  time.numb <- as.numeric(strsplit(time, "ZT")[[1]][[2]])
  if (!is.na(jrep)){
    time.numb.withrep <- time.numb + 24 * (jrep - 1)
  } else {
    time.numb.withrep <- time.numb - 48  # if NA, it is Liver, which starts at ZT50, make it start at ZT02 instead
  }
  add.leading.zeros <- FALSE
  if (add.leading.zeros){
    # add leading zeros (OPTIONAL)
    time.numb.withrep <- sprintf("%02d", time.numb.withrep)
    time.new <- paste0("ZT", time.numb.withrep)
  } else {
    time.new <- time.numb.withrep
  }
  return(time.new)
}, times, replicates, USE.NAMES = FALSE)

dat.bytranscript <- data.frame(gene = rep(genes, length(sampnames)),
                       tpm = unlist(dat),
                       tissue = rep(tissues, each = length(genes)),
                       time = rep(times.new, each = length(genes)),
                       geno = rep(genos, each = length(genes)),
                       feeding = rep(feedings, each = length(genes)),
                       experiment = "RNASeq")

dat.long <- dat.bytranscript %>%
  group_by(gene, tissue, time, geno, feeding, experiment) %>%
  summarise(exprs.linear = sum(tpm)) %>%
  filter(exprs.linear > 0)
dat.long$exprs <- log2(dat.long$exprs.linear + eps)

dat.long$geno <- factor(as.character(dat.long$geno), levels = c("SV129", "BmalKO"))
# Plot example ------------------------------------------------------------

jgene <- "Arntl"
jgene <- "Dbp"
jgene <- "Nr1d1"
jgene <- "Gm129"
jgene <- "Ciart"
jgene <- "Npas2"
ggplot(subset(dat.long, gene == jgene), aes(x = time, colour = tissue, linetype = geno, y = exprs)) + 
  geom_point() + geom_line() + xlab("ZT") + ylab("log2 expression") + 
  theme_bw(24) + 
  facet_wrap(~geno) + theme(aspect.ratio = 1, legend.position = "bottom")

PlotGeneAcrossTissues(subset(dat.long, gene == "Arntl"), "Arntl", make.pretty=TRUE, jxlab="ZT", do.facet.wrap=FALSE)

dir.create("Robjs/liver_kidney_atger_nestle")
save(dat.long, file="Robjs/liver_kidney_atger_nestle/dat.long.Robj")
