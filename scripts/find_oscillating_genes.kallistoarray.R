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

GetAmpRelToMax <- function(dat, fits.mostrhythmic){
  tissue <- as.character(dat$tissue[1])
  amp.max <- fits.mostrhythmic[tissue, ]$amp 
  dat$relamp <- dat$amp / amp.max
  return(dat)
}

FixColname <- function(tiss, time, exper){
  # Make MARA compatible colnames e.g., Adr22_array
  fixed_cname <- paste0(tiss, time, "_", exper)
  return(fixed_cname)
}

FixColnamesForMARA <- function(m){
  # Make colnames compatible with MARA
  # gene colname is "Gene.ID" (first column)
  # samples are TissueTime_arrayORrnaseq (e.g., Adr22_array)
  # 
  # Rearrrange colnames also to make put all arrays, THEN all rna-seq
  cnames <- colnames(m)
  genestr <- "Gene.ID"
  
  # Fix sample names (ignore first column which is gene)
  tissues <- sapply(cnames[2:length(cnames)], function(s) strsplit(s, "_")[[1]][[1]])
  times <- sapply(cnames[2:length(cnames)], function(s) strsplit(s, "_")[[1]][[2]])
  experiment <- sapply(cnames[2:length(cnames)], function(s) strsplit(s, "_")[[1]][[3]])
  
  sampnames <- mapply(FixColname, tissues, times, experiment)
  
  cnames.new <- c(genestr, sampnames)
  colnames(m) <- cnames.new
  
  # reorder colnames now
  array.samps <- cnames.new[grep("_array", cnames.new)]
  array.rnaseq <- cnames.new[grep("_rnaseq", cnames.new)]
  
  m <- m[, c(genestr, array.samps, array.rnaseq)]
  return(m)
}

WriteGeneListMat <- function(fits.rhyth.tiss, ka.long, outdir){
  # input: list of rhythmic genes by tissues (gene and tissue colnames), should be only ONE tissue (for dplyr)
  # dataframe of long matrix: (colnames: gene tissue time exprs experiment)
  # 
  # NOTE: expression should be centered across all samples! Otherwise it's a weird MARA model
  tiss <- unique(as.character(fits.rhyth.tiss$tissue))
  if (length(tiss) != 1){
    print('Warning: expect tissue length = 1')
  }
  fname.genelist <- paste(tiss, "genelist", sep = ".")
  fname.mat <- paste(tiss, "mat", sep=".")
  
  genes <- as.character(fits.rhyth.tiss$gene) 
  # write genelist
  sink(file = file.path(outdir, fname.genelist))
  for (g in genes){
    cat(g)
    cat("\n")
  }
  sink()
  
  # write matrix expression
  ka.sub <- subset(ka.long, gene %in% genes & tissue == tiss)
  m <- dcast(data = ka.sub, formula = gene ~ tissue + time + experiment, value.var = "exprs.centered")
  
  m <- FixColnamesForMARA(m)
  
  write.table(m, file = file.path(outdir, fname.mat), quote = FALSE, sep = '\t', col.names = TRUE, row.names = FALSE)
  return(data.frame(NULL))
}


# Load matrix -------------------------------------------------------------

array.path <- "data/exprs_matrices/array_adj_to_kallisto.slope07.txt"
rna.seq.path <- "/home/yeung/projects/tissue-specificity/data/kallisto/abundance.genecounts.matrix.txt"
ka.long <- LoadLong(array.path, rna.seq.path, scale.factor = 100, pseudocount = 1)

PlotGeneAcrossTissues(subset(ka.long, gene == "Arntl"))


# Center exprs for each gene across all samples ---------------------------

ka.centered <- ka.long %>%
  group_by(gene) %>%
  mutate(exprs.centered = scale(exprs, center = TRUE, scale = FALSE)

# Fit rhythmic ------------------------------------------------------------

fits <- FitRhythmicDatLong(ka.long)

save(fits, file = "Robjs/kallistoarray.fits.Robj")


# P-value adjuist ---------------------------------------------------------

fits$pval.adj <- NULL
fits <- fits %>%
  group_by(tissue) %>%
  mutate(pval.adj = p.adjust(pval)) %>%
  arrange(pval.adj)
# fits$pval.adj = p.adjust(fits$pval)

# Find most rhythmic gene: use this as relative amplitude for othe --------

max.pval <- 1e-3
# max.pvaladj <- 0.9

fits.mostrhythmic <- fits %>%
  group_by(tissue) %>%
  filter(pval <= max.pval) %>%
#   filter(pval.adj <= max.pvaladj) %>%
  do(FindMostRhythmic(.)) %>%
  data.frame(.)

rownames(fits.mostrhythmic) <- fits.sort$tissue  # indexing


# Get amplitude relative to max amp ---------------------------------------

fits.relamp <- fits %>%
  group_by(tissue) %>%
  do(GetAmpRelToMax(., fits.mostrhythmic))

# List rhythmic genes -----------------------------------------------------

min.relamp <- 0.1  # fraction of largest amplitude of most rhythmic gene

fits.rhyth <- subset(fits.relamp, pval <= max.pval & amp >= min.relamp)
# fits.rhyth <- subset(fits.relamp, pval <= max.pval & relamp >= min.relamp)
# fits.rhyth <- subset(fits.relamp, pval.adj <= max.pvaladj & relamp >= min.relamp)

rhythmic.genes <- unique(fits.rhyth$gene)

print(length(rhythmic.genes))

# number of rhythmic genes in each tissue
rhyth.by.tiss <- fits.rhyth %>%
  group_by(tissue) %>%
  summarise(length(gene))

print(rhyth.by.tiss)


# Print gene list and expression matrix -----------------------------------

outdir <- paste0("/home/yeung/projects/tissue-specificity/data/gene_lists/rhythmic_genes_by_tissues_kallistoarray", ".pval", max.pval, ".relamp", min.relamp)
dir.create(outdir)

fits.rhyth %>%
  group_by(tissue) %>%
  do(WriteGeneListMat(., ka.centered, outdir))
