# 2015-06-18
# find_oscillating_genes.kallistoarray.R
# We want to find oscillating genes, with p-value and amplitude proportional 
# to a core clock gene as a cutoff.

library(ggplot2)
library(dplyr)
# Functions ---------------------------------------------------------------

source("scripts/functions/PlotGeneTpm.R")
source("scripts/functions/LoadArray.R")
source("scripts/functions/LoadKallistoGene.R")
source("scripts/functions/GetTissueTimes.R")
source("scripts/functions/FitRhythmic.R")
source("scripts/functions/PlotGeneAcrossTissues.R")

FindMostRhythmic <- function(dat, colname="amp", decreasing = TRUE){
  # return first row after sorting by colname
  return(dat[order(dat[[colname]], decreasing = decreasing), ][1, ])
}

LoadLong <- function(array.path, rna.seq.path, scale.factor = 1, pseudocount = 1e-5){
  kallisto.wide <- LoadKallistoGene(rna.seq.path, form = "wide")  # adjusts colnames to match array
  array.wide <- read.table(array.path)
  # remove rows with negatives (should be 7 of them)
  problem.genes <- rownames(array.wide[which(apply(array.wide, 1, min) < 0), ])
  print("Problem genes:")
  print(problem.genes)
  good.genes <- setdiff(rownames(array.wide), problem.genes)
  array.wide <- array.wide[good.genes, ]
  common.genes <- intersect(rownames(kallisto.wide), rownames(array.wide))
  kallisto.sub <- kallisto.wide[common.genes, ]
  array.sub <- array.wide[common.genes, ]
  
  tissues.rnaseq <- GetTissues(colnames(kallisto.sub), get_unique = FALSE)
  times.rnaseq <- GetTimes(colnames(kallisto.sub), get_unique = FALSE)
  tissues.array <- GetTissues(colnames(array.sub), get_unique = FALSE)
  times.array <- GetTimes(colnames(array.sub), get_unique = FALSE)
  
  ka.long <- data.frame(gene = c(rep(rownames(kallisto.sub), ncol(kallisto.sub)), rep(rownames(array.sub), ncol(array.sub))),
                        tissue = c(rep(tissues.rnaseq, each = nrow(kallisto.sub)), rep(tissues.array, each = nrow(array.sub))),
                        time = as.numeric(c(rep(times.rnaseq, each = nrow(kallisto.sub)), rep(times.array, each = nrow(array.sub)))),
                        exprs = c(unlist(kallisto.sub), unlist(array.sub)),
                        experiment = c(rep("rnaseq", nrow(kallisto.sub) * ncol(kallisto.sub)), rep("array", nrow(array.sub) * ncol(array.sub))))
  ka.long$exprs <- log2(scale.factor * ka.long$exprs + pseudocount)
  return(ka.long)
}

PlotGeneAcrossTissues2 <- function(dat, jgene){
  dat.sub <- subset(dat, gene == jgene)
  m <- ggplot(dat.sub, aes(x = time, y = tpm, colour = experiment, group = experiment, fill = experiment)) +
    geom_point() + geom_line() + facet_wrap(~tissue)
    ggtitle(jgene)
  print(m)
}

# Load matrix -------------------------------------------------------------

array.path <- "data/exprs_matrices/array_adj_to_kallisto.slope07.txt"
rna.seq.path <- "/home/yeung/projects/tissue-specificity/data/kallisto/abundance.genecounts.matrix.txt"
ka.long <- LoadLong(array.path, rna.seq.path, scale.factor = 100, pseudocount = 1)

PlotGeneAcrossTissues(subset(ka.long, gene == "Arntl"))


# Fit rhythmic ------------------------------------------------------------

fits <- FitRhythmicDatLong(ka.long)


# Find most rhythmic gene: use this as relative amplitude for othe --------

max.pval <- 1e-5
fits.mostrhythmic <- fits %>%
  group_by(tissue) %>%
  filter(pval <= max.pval) %>%
  do(FindMostRhythmic(.)) %>%
  data.frame(.)

rownames(fits.mostrhythmic) <- fits.sort$tissue  # indexing


# Get amplitude relative to max amp ---------------------------------------

GetAmpRelToMax <- function(dat, fits.mostrhythmic){
  tissue <- as.character(dat$tissue[1])
  amp.max <- fits.mostrhythmic[tissue, ]$amp 
  dat$relamp <- dat$amp / amp.max
  return(dat)
}

fits.relamp <- fits %>%
  group_by(tissue) %>%
  do(GetAmpRelToMax(., fits.mostrhythmic))

# List rhythmic genes -----------------------------------------------------

min.relamp <- 0.3  # fraction of largest amplitude of most rhythmic gene

fits.rhyth <- subset(fits.relamp, pval <= max.pval & relamp >= min.relamp)

rhythmic.genes <- unique(fits.rhyth$gene)

print(length(rhythmic.genes))

# number of rhythmic genes in each tissue
rhyth.by.tiss <- fits.rhyth %>%
  group_by(tissue) %>%
  summarise(length(gene))

print(rhyth.by.tiss)
