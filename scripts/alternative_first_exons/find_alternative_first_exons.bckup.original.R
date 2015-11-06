# 2015-04-28
# find_alternative_first_exons.R

setwd("~/projects/tissue-specificity/")

library(ggplot2)
library(plyr)
library(mixtools)
library(doParallel)
# Functions ---------------------------------------------------------------

source("scripts/functions/LoadArrayRnaSeq.R")
source("scripts/functions/LoadAndHandleData.R")
source("scripts/functions/MixtureModelFunctions.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/MakeCluster.R")
source("scripts/functions/ReadListToVector.R")
source("scripts/functions/GrepRikGenes.R")

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

ShannonEntropy <- function(x.vec, normalizeout=TRUE){
  # should sum to 1
  x.vec <- x.vec / sum(x.vec)
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

# dat.long <- LoadArrayRnaSeq()
dat <- LoadRnaSeq()


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

# cov.long.filt$reads_norm <- cov.long.filt$reads / cov.long.filt$rnaseq_reads
cov.long.split <- split(cov.long, cov.long$gene)


# Average across tissues --------------------------------------------------

# rscript_args <- MakeCluster()
# ncores <- detectCores()
# cl <- makeCluster(ncores, rscript_args=rscript_args)
# # stopCluster(cl)  # if you cannot open connections
# registerDoParallel(cl)
doParallel::registerDoParallel(cores = 48)
start <- Sys.time()
cov.long.filt.avg <- ddply(cov.long.filt, 
                           .(tissue, transcript),
                           .fun = summarise,
                           reads_norm_avg = mean(reads_norm),
                           gene = unique(gene),
                           .parallel = TRUE)
print(Sys.time() - start)

# save(cov.long.filt.avg, file = "data/alternative_exon_usage/cov.long.filt.avg.Robj")
# load("data/alternative_exon_usage/cov.long.filt.avg.Robj")

# Calculate ShannonEntropy ------------------------------------------------

cov.long.filt.entropy <- ddply(cov.long.filt.avg, .(transcript),
                               .fun = summarise,
                               entropy = ShannonEntropy(reads_norm_avg))


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
