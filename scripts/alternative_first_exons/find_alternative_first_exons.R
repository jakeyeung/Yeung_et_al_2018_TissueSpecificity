# 2015-04-28
# find_alternative_first_exons.R

setwd("~/projects/tissue-specificity/")

library(ggplot2)
library(dplyr)
library(mixtools)
library(gplots)
library(parallel)
# Functions ---------------------------------------------------------------

source("scripts/functions/LoadArrayRnaSeq.R")
source("scripts/functions/LoadAndHandleData.R")
source("scripts/functions/MixtureModelFunctions.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/MakeCluster.R")
source("scripts/functions/ReadListToVector.R")
source("scripts/functions/GrepRikGenes.R")
source("scripts/functions/FitRhythmic.R")

cossim <- function(x, y){
  return(x %*% y / sqrt(x%*%x * y%*%y))
}

LoopCor <- function(m, show.which=FALSE, input.vec1=NA){
  imin <- NA
  jmin <- NA
  jcor.min <- 2  # init because pearson cor is between -1 and 1. Use a number outside of this range.
  for (i in 1:ncol(m)){
    if (!is.na(input.vec1)){
      i <- input.vec1
    }
    vec1 <- m[, i]
    for (j in i:ncol(m)){
      if (j == i){
        next  # dont need to compare between same vec
      }
      vec2 <- m[, j]
      jcor <- cossim(vec1, vec2)
      if (jcor < jcor.min){
        jcor.min <- jcor
        imin <- i
        jmin <- j
      }
    }
    if (!is.na(input.vec1)){
      break
    }
  }
  if (show.which){
    jtissues <- c(colnames(m)[imin], colnames(m)[jmin])
    print(jtissues)
    return(list(tissues = jtissues,
                min.cor = jcor.min))
  } 
  return(jcor.min)
}

GetMinCor <- function(df){
  m <- acast(data = df, transcript ~ tissue, value.var = "norm_reads")
  jcor.min <- LoopCor(m)
  return(data.frame(min.cor = jcor.min))
}

Normalize <- function(x, pseudocount = 1){
  if (pseudocount > 0){
    x <- x + pseudocount
  }
  x.norm <- x / sum(x)
  return(x.norm)
}

AvgAcrossTranscripts <- function(df){
  ddply(df, .(transcript, gene, tissue), summarise, mean_reads = mean(reads))
}

NormalizeReads <- function(df){
  # Normalize reads across transcripts
  ddply(df, .(gene), transform, norm_reads = Normalize(mean_reads))
}

GetGeneNameFromAnnot <- function(annot){
  # from gene_name=RP23-271O17.1;transcript_id=ENSMUST00000193812 retrieve RP23-27017.1
  gene.str <- strsplit(as.character(annot), ';')[[1]][[1]]
  gene.str <- strsplit(gene.str, '=')[[1]][[2]]
  return(gene.str)
}

GetTranscriptIDFromAnnot <- function(annot){
  # from gene_name=RP23-271O17.1;transcript_id=ENSMUST00000193812 retrieve ENSMUST...
  transcript.str <- strsplit(as.character(annot), ';')[[1]][[2]]
  transcript.str <- strsplit(transcript.str, '=')[[1]][[2]]
  return(transcript.str)
}

rowMax <- function(df){
  # Return vector of maximums from df
  return(apply(df, 1, max))
}

GetTissuesAFE <- function(x){
  # Get tissues from column names from Adr_CT22 format
  substr(x, "_")[[1]][[1]]
}

TissueMapping <- function(cov.to.rnaseq = TRUE){
  # Tissue names from coverage are slightly different from
  # tissue names from RNASeq data. Create the mapping between 
  # coverage to rnaseq tissue names
  list("Adr" = "Adr",
       "Aor" = "Aorta",
       "BFat" = "BFAT",
       "Bstm" = "BS",
       "Cer" = "Cere",
       "Hrt" = "Heart",
       "Hyp" = "Hypo",
       "Kid" = "Kidney",
       "Liv" = "Liver",
       "Lun" = "Lung",
       "Mus" = "Mus",
       "WFat" = "WFAT")
}

ConvertTissuesToMatchRnaSeq <- function(tissues){
  # make tissue names look like RNASeq column names using TissueMapping
  tissue.map <- TissueMapping()
  sapply(tissues, function(x){
    tissue.map[[x]]
  })
}

GetExprsAsVector <- function(dat, genes, tissuetime){
  # Given a vector of rownames and column names, extract
  # its corresponding element in the dat. Return as
  # a vector.
  if (length(genes) != length(tissuetime)) print("Genes and tissuetime is not same length")
  lookups <- mapply(function(x, y){
    return(dat[x, y])
  }, x = genes, y = tissuetime)
  return(lookups)
}

GetLocationFromAnnotation <- function(bed, gene_name, transcript_id){
  # Given bed, gene_name and transcript_id, return the chromo, start, end
  annot <- paste0("gene_name=", gene_name, ";transcript_id=", transcript_id)
  sub <- subset(bed, annotations == annot)
  # return as UCSC-style
  return(paste0(sub$chromosome, ":", sub$start, "-", sub$end))
}

SubsetBed <- function(bed, gene_name, transcript){
  # Subset bed based on grepping annotations from gene name
  if (missing(transcript)){
    return(bed[grepl(gene_name, bed$annotations), ])  
  } else if (missing(gene_name)){  
    return(bed[grepl(transcript, bed$annotations), ])  
  }
}

NormalizeBySum <- function(x){
  # Normalize a vector by its sum
  return(x / sum(x))
}

GetFirst <- function(x){
  return(x[1])
}

ShannonEntropy <- function(x.vec, normalize=FALSE){
  if (normalize){
      # should sum to 1
    x.vec <- x.vec / sum(x.vec)
  }
  entropy <- 0
  for (x in x.vec){
    entropy <- entropy + x * log2(1 / x)
  }
  entropy <- entropy / log2(length(x.vec))
  return(entropy)
}

MulitpleStarts <- function(df, min_dist = 500){
  # Check if df has multiple exons, ddply from bed
  if (nrow(df) <= 1){
    return(data.frame(MultiStart = FALSE))
  }
  dist <- diff(c(min(df$start), max(df$start)))
  if (dist >= min_dist){
    return(data.frame(MultiStart = TRUE))
  } else{
    return(data.frame(MultiStart = FALSE))
  }
}

# Load RNASeq ---------------------------------------------------

dat.long <- LoadArrayRnaSeq()
dat <- LoadRnaSeq()


# Find rhythmic genes -----------------------------------------------------

dat.long.by_genetiss <- group_by(dat.long, gene, tissue)
dat.long.by_genetiss.split <- split(dat.long.by_genetiss, dat.long.by_genetiss$tissue)
print("Finding rhythmic genes (~3 minutes)")
start <- Sys.time()
dat.fitrhyth.split <- mclapply(dat.long.by_genetiss.split, function(jdf){
  do(.data = jdf, FitRhythmic(df = .))
}, mc.cores = 12)
dat.fitrhyth <- do.call(rbind, dat.fitrhyth.split)
print(Sys.time() - start)
if (exists("dat.fitrhyth")) rm(dat.fitrhyth.split)  # GC

# Find cutoff for expressed genes -----------------------------------------

# exprs.vec <- log2(unlist(dat)[which(unlist(dat) > 0)])
# plot(density(exprs.vec))
# 
# # mixture of two gaussians?
# mixmdl <- normalmixEM(exprs.vec, lambda = c(0.25, 0.75), mu = c(2.5, 9), k = 2)
# plot(mixmdl,which=2)
# lines(density(exprs.vec), lty=2, lwd=2)
# 
# cutoff <- optimize(ShannonEntropyMixMdl, interval = c(2, 8), mixmdl = mixmdl, maximum = TRUE)
# cutoff <- cutoff$maximum  # cutoff = 4.883356
cutoff <- 4.883356

print(paste("Cutoff:", cutoff))

# Main --------------------------------------------------------------------


# Load data fix rownames ---------------------------------------------------------------


infile <- "data/alternative_exon_usage//Mus_musculus.GRCm38.79.first_exons.coveragebed"

exon.cov <- read.table(infile, header = TRUE, sep = "\t")

gene.names <- sapply(exon.cov$annotations, GetGeneNameFromAnnot)
transcript.ids <- sapply(exon.cov$annotations, GetTranscriptIDFromAnnot)

rnames <- make.names(gene.names, unique = TRUE)

rownames(exon.cov) <- rnames


# Define annotations and data matrix --------------------------------------


# Make a matrix of only coverage values
bed <- exon.cov[, 1:6]  # we might need it later
exon.cov <- exon.cov[, 7:ncol(exon.cov)]

plot(density(log2(unlist(exon.cov))))




# Match genes with RNA-Seq ------------------------------------------------

dat.filt <- dat[make.names(gene.names), ]  # handle these X5830428H23Rik-like gnes
rownames(dat.filt) <- rownames(exon.cov)

# Convert to long ---------------------------------------------------------

# Get Adr from Adr_CT22
tissues <- sapply(colnames(exon.cov), function(x) strsplit(x, "_")[[1]][[1]])
# Convert mappings to match RNASeq tissuenames
tissues <- ConvertTissuesToMatchRnaSeq(tissues)
# Get 22 from Adr_CT22
times <- sapply(colnames(exon.cov), function(x) strsplit(x, "_")[[1]][[2]])
times <- as.numeric(sapply(times, function(x) substr(x, nchar(x) - 1, nchar(x))))

cov.long <- data.frame(gene = rep(gene.names, ncol(exon.cov)),
                       transcript = rep(transcript.ids, ncol(exon.cov)),
                       tissue = rep(tissues, each = length(gene.names)),
                       time = rep(times, each = length(gene.names)),
                       reads = unlist(exon.cov),
                       rnaseq_reads = unlist(dat.filt))
head(cov.long)


# Only consider tissue-specific rhythmic genes ----------------------------

# from get_rhythmic_genes_from_list.R
rhythmic_genes <- ReadListToVector("results/tissue_specific_rhythmic_genes/tissue_specific_rhythmic_genes.genome_wide.cut.txt", 
                                   HEADER = TRUE)
# Remove "X" from Rik genes
Xgenes <- GrepRikGenes(gene.list = rhythmic_genes)
Xgenes <- RemoveX(Xgenes) 

# Find transcripts whose start sites are significantly far away -----------

# decrease number of genes to consider
common.genes <- intersect(unique(gene.names), c(rhythmic_genes, Xgenes))  # 1062 genes

bed$gene <- gene.names
bed.filtered <- ddply(bed, .(gene), MulitpleStarts)

gene.multistarts <- bed.filtered$gene[which(bed.filtered$MultiStart == TRUE)]

# Furhter filter genes
common.genes <- intersect(common.genes, gene.multistarts)


# Filter rhythmic genes with multistarts ----------------------------------

# cov.long.filts <- subset(cov.long, gene %in% common.genes)

# Filter lowly expressed genes --------------------------------------------

# cov.long.filt <- subset(cov.long.filts, rnaseq_reads >= 2^cutoff)  # cutoff established in log2 scale


# Optionaly do not filter anything ----------------------------------------

cov.long.filt <- cov.long


# Normalize exon reads ----------------------------------------------------

# cov.long.filt$reads_norm <- cov.long.filt$reads / cov.long.filt$rnaseq_reads  # naive

# cov.normreads <- cov.long %>%
#   group_by(tissue, gene, time) %>%
#   mutate(norm_reads = Normalize(reads), n_starts = length(reads))

by_tissuegene <- group_by(cov.long, transcript, tissue, gene)
cov.avgreads <- summarise(by_tissuegene, mean_reads = mean(reads))
cov.avgreads.by_tissuegene <- group_by(cov.avgreads, tissue, gene)
cov.normreads <- mutate(cov.avgreads.by_tissuegene, norm_reads = Normalize(mean_reads), n_starts = length(mean_reads))


# Remove n_starts == 1 ----------------------------------------------------

cov.normreads <- subset(cov.normreads, n_starts > 1)


# Remove lowly expressed genes --------------------------------------------

# find cutoff
normreads.vec <- log2(cov.normreads$mean_reads + 1)
mixmdl.normreads <- normalmixEM(normreads.vec, lambda = c(0.5, 0.5), mu = c(0.1, 6), k = 2)
plot(mixmdl.normreads, which = 2)
lines(density(normreads.vec), lty = 2, lwd = 2)
cutoff.normreads <- optimize(ShannonEntropyMixMdl, interval = c(1, 5), mixmdl = mixmdl.normreads, maximum = TRUE)
cutoff.normreads <- 2^(cutoff.normreads$maximum)
print(paste("mean reads cutoff:", cutoff.normreads))  # 1.01

cov.normreads.filt <- mutate(cov.normreads, mean_reads.gene = mean(mean_reads))

cov.normreads.filt <- subset(cov.normreads.filt, mean_reads.gene > cutoff.normreads)

# Calculate maximum difference --------------------------------------------

cov.normreads.sub <- subset(cov.normreads.filt, gene %in% common.genes)
cov.normreads.by_gene <- group_by(cov.normreads.sub, gene)

# find AFEs by "minimum correlation"
cov.mincor <- do(.data = cov.normreads.by_gene, GetMinCor(df = .))  # super slow


# Calculate log2 fold change ----------------------------------------------

# find AFEs by log2 fold change between "rhythmic" and "not rhythmic" genes
cov.normreads.by_genetiss <- group_by(cov.normreads, gene, tissue)
dat.fitrhyth.filt <- data.frame(subset(dat.fitrhyth, gene %in% common.genes, select = c(gene, tissue, as.numeric(pval), amp)))  # faster
rnames <- apply(dat.fitrhyth.filt, 1, function(x) paste0(x[2], '-', x[1]))
rownames(dat.fitrhyth.filt) <- rnames; rm(rnames)
dat.fitrhyth.filt <- subset(dat.fitrhyth.filt, select = c(pval, amp))
dat.fitrhyth.filt <- data.matrix(dat.fitrhyth.filt)  # faster to work with matrices

# ask if rhythmic 20 seconds
rhythmic.or.not <- apply(cov.normreads.by_gene, 1, GetRhythmicOrNot, fitdf = dat.fitrhyth.filt)

rhythmic.or.not <- cov.isrhyth.vec

# Append to df
cov.normreads.by_gene.rhyth <- cbind(cov.normreads.by_gene, rhythmic.or.not)

# Find 

GetRhythmicOrNot <- function(x, fitdf){
  # Expect x to be a row from cov.normreads, with tissue and gene in 2nd and 3rd col
  tiss <- x[2]
  gene <- x[3]
  rname <- paste0(tiss, '-', gene)
  fitdf.sub = tryCatch({
    fitdf[rname, ]
  }, warning = function(w) {
    print("Warning")
    print(w)
  }, error = function(e) {
    # print(paste("Cannot access:", rname))
    return(NA)
  })  
  if (is.na(fitdf.sub[1])){
    return(NA)
  }
  pval <- fitdf.sub[1]
  amp <- fitdf.sub[2]
  annots <- RhythmicOrNot(pval, amp)
  return(annots)
}

RhythmicOrNot <- function(pval, amp, min.pval = 1e-5, max.pval = 0.05, max.amp = 0.5, min.amp = 0.1){
  if (pval < min.pval & amp > max.amp){
    return("Rhythmic")
  } else if (pval > max.pval & amp < min.amp){
    return("NotRhythmic")
  } else {
    return(NA)  # undecided
  }
}

FoldChangeRhyth <- function(jdf){
  # Calculate log2 fold change between "rhythmic" and "non rhythmic" genes 
}



# Show distributions ------------------------------------------------------

(head(as.data.frame(cov.mincor[order(cov.mincor$min.cor), ]), n = 50))
plot(density(cov.mincor$min.cor), xlim=c(0, 1))


# Do sanity checks on examples --------------------------------------------

jgene <- "Ddc"
jgene <- "Sgk2"
jgene <- "Gas7"
jgene <- "Insig2"
# plot exprs
PlotGeneAcrossTissues(subset(dat.long, gene == jgene))
# check normreads
jdf <- subset(cov.normreads.filt, gene == jgene)
print(data.frame(jdf))
# check out the matrix for correlations
m <- acast(jdf, transcript ~ tissue, value.var = "norm_reads")
# check which two tissues were called as "most different"
best.cor <- LoopCor(m, show.which = TRUE)
# check which transcript accounts for largest difference
(m.sub <- m[, best.cor$tissues])
# return transcript with highest difference
m.diffs <- apply(m.sub, 1, function(x) abs(log2(x[1] / x[2])))
(m.diffs.max <- m.diffs[which(m.diffs == max(m.diffs))])
transcript.max <- names(m.diffs.max)
# get chromosome location of this transcript
GetLocationFromAnnotation(bed, gene_name = jgene, transcript_id = transcript.max)
# Heatmap the matrix of correlations
my_palette <- colorRampPalette(c("white", "black"))(n = 299)
heatmap.2(m, density.info = "density", trace = "none", margins = c(8, 14), main = jgene, col = my_palette, cexRow=1.1)

MaxDiff <- 

# Calculate ShannonEntropy ------------------------------------------------

cov.entropy <- summarise(cov.normreads, entropy = ShannonEntropy(norm_reads))

# Remove NaNs -------------------------------------------------------------

cov.long.filt.entropy <- cov.long.filt.entropy[which(!is.nan(cov.long.filt.entropy$entropy)), ]

# Sort by entropy ---------------------------------------------------------

print(head(cov.long.filt.entropy[order(cov.long.filt.entropy$entropy), ], n = 50))

# Test with Ddc and Insig2 ------------------------------------------------

jgene <- "Dbp"
jtrans <- "ENSMUST00000107740"

jgene <- "Insig2"
jtrans <- "ENSMUST00000161068"

jgene <- "Ddc"
jtrans <- "ENSMUST00000178704"

jgene <- "Hnf4a"
jtrans <- "ENSMUST00000143911"

cov.long.sub <- subset(cov.long, gene == jgene)

test <- subset(cov.long.sub, transcript == jtrans)
test$normcov <- test$reads / test$rnaseq_reads

ggplot(data = test, aes(x = time, y = reads)) + 
  geom_point() + 
  geom_line() + 
  facet_wrap(~tissue) + 
  ggtitle(paste(jgene, jtrans))
