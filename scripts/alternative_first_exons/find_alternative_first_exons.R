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

SubsetBed <- function(bed, gene_name){
  # Subset bed based on grepping annotations from gene name
  return(bed[grepl(gene_name, bed$annotations), ])
}

NormalizeBySum <- function(x){
  # Normalize a vector by its sum
  return(x / sum(x))
}

GetFirst <- function(x){
  return(x[1])
}

ShannonEntropy <- function(x.vec, scaleit=TRUE, normalizeout=TRUE){
  if (scaleit){
    x.vec <- x.vec / sum(x.vec)
  }
  entropy <- 0
  for (x in x.vec){
    if (x > 0){
      entropy <- entropy + x * log2(1 / x)
    }
  }
  if (normalizeout){
    # Divide by max, log2(n)
    entropy <- entropy / log2(length(x.vec))
  }
  return(entropy)
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

infile <- "data/alternative_exon_usage//Mus_musculus.GRCm38.79.first_exons.coveragebed"

exon.cov <- read.table(infile, header = TRUE, sep = "\t")

gene.names <- sapply(exon.cov$annotations, GetGeneNameFromAnnot)
transcript.ids <- sapply(exon.cov$annotations, GetTranscriptIDFromAnnot)

rnames <- make.names(gene.names, unique = TRUE)

rownames(exon.cov) <- rnames

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



# Filter lowly expressed genes --------------------------------------------

cov.long.filt <- subset(cov.long, rnaseq_reads >= 2^cutoff)  # cutoff established in log2 scale


# Normalize exon reads ----------------------------------------------------

cov.long.filt$reads_norm <- cov.long.filt$reads / cov.long.filt$rnaseq_reads


# Average across tissues --------------------------------------------------
# ** Ddc: clear alternative promoter usage Kidney and Liver rhythmic but not in others.
# Adam19: example of differential promoter usage (the promoters are close together). Rhythmic in Heart
# Ahcyl2: a small alternative promoter is used in lung. Not sure if it's definitive. Rhythmic in Lung
# Chka: no significant alternative promoter. RHythmic in Liver.
# Slc41a3: alternative promoter used in Adr and Heart. But only heart is rhythmic. Example of both alternative and differential promoter usage?
# Slc6a6: Rhythmic in WFAT, BFAT, LUNG, Aorta, but no sign of alterntaive promoter usage. Seems it was annotated incorrectly, the second promoter doesn't have any annotated start site.
# Angptl2: misannotated? Rhythmic in many tissues but no sign of alternative promoter usage.
# Tshr: no alt promoter usage
# Pxmp4: no alt promoter usage
# Pnp: no alt promoter usage
# Sec14l2: curious Adr triple peak signal, all others are flat. But UCSC doesnt show alt promoter
# Hdac5: rhythmic in liver only but UCSC doesnt show alt promoter. Unclear.
# Dtna: "rhyhmic in BFAT, but no sign that this is due to alt promoter. Alt promoter usage noticeable in BS, Cere and Hypo and Lung, however. Differences between heart and BS is obvious.
# Ppargc1a: small evidence of alt promoter usage. 
# Ank3: kidney and lung small evidence of alt promoter usage but doesn't match expression profile.
# ** Sept9: some evidence that BFAT and Liver take the upstream promoter, which are the rhythmic genes.
# Ncoa7: small evidence of alternative promoter usage linked with rhythmicity. But rhythms are small.
# Asph: strong evidence of alternative promoter usage (see it in muscle) but BFAT has nothing particularly special.
# Pnkd: sign that liver may be alternative promoter usage, but it's weak
# Ngef: sign of alternative promoter usage for kidney and liver, but kidney and liver is weakly rhythmic. Lung is strongly rhythmic.
# ** Steap3: Liver specific promoter.
# ** Insig2 clear alt promoter in liver. Strong rhythmics in liver.
# ** Slc45a3 strong rhythmic in liver and alt promoter in liver.
# ** Mpzl1: strong rhythmic in liver, small rhythms in other genes. Possible alt promoter? Should check.

# genes.test <- c("Arntl", "Clock", "Insig2", "Dbp", "Ddc", "Sept9", "Steap3", "Slc45a3", "Mpzl1", "Dtna", "Tshr", "Pnp",
#                 "Hdac5", "Ppargc1a", "Ank3", "Ncoa7", "Hnf4a", "Nr1d1", "Rest")
# rscript_args <- MakeCluster()
# ncores <- detectCores()
# cl <- makeCluster(ncores, rscript_args=rscript_args)
# # stopCluster(cl)  # if you cannot open connections
# registerDoParallel(cl)
doParallel::registerDoParallel(cores = 40)
start <- Sys.time()
cov.long.filt.avg <- ddply(cov.long.filt, 
                           .(tissue, transcript),
                           .fun = summarise,
                           reads_norm_avg = mean(reads_norm),
                           gene = GetFirst(gene),
                           .parallel = TRUE)
print(Sys.time() - start)

save(cov.long.filt.avg, file = "data/alternative_exon_usage/cov.long.filt.avg.Robj")

# # Force their sums to 1 ---------------------------------------------------
# 
# cov.long.filt.avg$frac_reads_norm <- NormalizeBySum(cov.long.filt.avg$reads_norm_avg)
# 
# 
# # Calculate ShannonEntropy ------------------------------------------------
# 
# cov.long.filt.entropy <- ddply(cov.long.filt.avg, .(transcript),
#                                summarise,
#                                entropy = ShannonEntropy(frac_reads_norm, normalize = TRUE),
#                                gene = GetFirst(gene))
# 
# 
# # Sort by entropy ---------------------------------------------------------
# 
# print(cov.long.filt.entropy[order(cov.long.filt.entropy$entropy), ])
# 
# # Test with Ddc and Insig2 ------------------------------------------------
# 
# jgene <- "Dbp"
# jtrans <- "ENSMUST00000107740"
# 
# jgene <- "Insig2"
# jtrans <- "ENSMUST00000161068"
# 
# jgene <- "Ddc"
# jtrans <- "ENSMUST00000178704"
# 
# jgene <- "Hnf4a"
# jtrans <- "ENSMUST00000143911"
# 
# cov.long.sub <- subset(cov.long, gene == jgene)
# 
# test <- subset(cov.long.sub, transcript == jtrans)
# test$normcov <- test$reads / test$rnaseq_reads
# 
# ggplot(data = test, aes(x = time, y = reads)) + 
#   geom_point() + 
#   geom_line() + 
#   facet_wrap(~tissue) + 
#   ggtitle(paste(jgene, jtrans))
