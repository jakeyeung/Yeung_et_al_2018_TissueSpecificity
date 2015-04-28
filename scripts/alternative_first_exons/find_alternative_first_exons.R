# 2015-04-28
# find_alternative_first_exons.R

setwd("~/projects/tissue-specificity/")

# Functions ---------------------------------------------------------------

source("scripts/functions/LoadAndHandleData.R")

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

# Load RNASeq and Array ---------------------------------------------------

dat <- LoadRnaSeq()

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
                       transcript.ids = rep(transcript.ids, ncol(exon.cov)),
                       tissues = rep(tissues, each = length(gene.names)),
                       times = rep(times, each = length(gene.names)),
                       reads = unlist(exon.cov),
                       rnaseq = unlist(dat.filt))
head(cov.long)


# # Get RNA-Seq gene expression ---------------------------------------------
# 
# # DESEQ2-normalized counts
# genes <- cov.long$gene
# tissuetime <- paste0(cov.long$tissues, cov.long$times)
# # takes about 15 minutes I think
# start <- Sys.time()
# # frac <- 0.05
# # exprs.rnaseq <- GetExprsAsVector(dat, genes[1:(frac * length(genes))], tissuetime[1:(frac * length(genes))])
# exprs.rnaseq <- GetExprsAsVector(dat, genes, tissuetime)
# print(Sys.time() - start)
# 
# save(exprs.rnaseq, file = "results/alternative_exon_usage/exprs.rnaseq.vector.Robj")
# cov.long$exprs <- exprs.rnaseq
# save(cov.long, file = "results/alternative_exon_usage/cov.long.Robj")
# print("DONE")

# Test with Ddc and Insig2 ------------------------------------------------

# jgene <- "Insig2"
# cov.long.sub <- subset(cov.long, gene == jgene)
# 
# subset(cov.long.sub, transcript.ids == "ENSMUST00000159528")
