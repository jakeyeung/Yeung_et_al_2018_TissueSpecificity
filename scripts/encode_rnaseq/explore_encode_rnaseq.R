# Jake Yeung
# Look at Kallisto RNA-Seq data from ENCODE
# 2015-08-12

library(wordcloud)
library(dplyr)
library(ggplot2)

source("scripts/functions/SortByTissue.R")

LoadEncodeRnaseq <- function(exprs.path, long.format = TRUE, transform.log2 = TRUE, scale.factor = 1000, pseudocount = 1){
  # Load data and put it into long format
  # if not long.format, then returns just the signal matrix
  # if log2, transforms by log2(tpm * scale.factor + pseudocount)
  
  exprs.dat <- read.table(exprs.path, header = TRUE)
  
  exprs.annot <- exprs.dat[, 1:6]
  exprs.signal <- exprs.dat[, 7:ncol(exprs.dat)]
  
  if (long.format == FALSE){
    if (transform.log2){
      exprs.signal <- log2(exprs.signal * scale.factor + pseudocount)
    }
    return(exprs.signal)
  }
  
  # capitalize Tissues because and remove underscores
  tissues <- unname(sapply(colnames(exprs.signal), 
                           function(sampid) CapitalizeRemoveUnderscores(strsplit(sampid, "\\.")[[1]][[1]])))
  
  bioreps <- unname(sapply(colnames(exprs.signal), 
                           function(sampid){
                             biorep_str <- strsplit(sampid, "\\.")[[1]][[2]]
                             biorep_i <- strsplit(biorep_str, "biorep")[[1]][[2]]
                           }))
  
  techreps <- unname(sapply(colnames(exprs.signal), 
                            function(sampid){
                              techrep_str <- strsplit(sampid, "\\.")[[1]][[3]]
                              techrep_i <- strsplit(techrep_str, "techrep")[[1]][[2]]
                            }))
  
  exprs.long <- data.frame(gene = rep(exprs.dat$gene_name, ncol(exprs.signal)),
                           transcript = rep(exprs.dat$target_id, ncol(exprs.signal)), 
                           tissue = rep(tissues, each = nrow(exprs.dat)),
                           biorep = rep(bioreps, each = nrow(exprs.dat)),
                           techrep = rep(techreps, each = nrow(exprs.dat)),
                           # sampleid = rep(colnames(exprs.signal), each = nrow(exprs.dat)),  # as sanity check
                           tpm = unlist(exprs.signal))
  
  # take mean of tech reps, then bio reps
  exprs.long.mean <- exprs.long %>%
    group_by(gene, transcript, tissue, biorep) %>%
    summarise(tpm = mean(tpm)) %>%
    group_by(gene, transcript, tissue) %>%
    summarise(tpm = mean(tpm))
  
  exprs.long.mean$logtpm <- log2(exprs.long.mean$tpm * scale.factor + pseudocount)
  
  return(exprs.long.mean)
}

Capitalize <- function(s){
  # capitalize first letter of string
  first.letter <- toupper(substr(s, 1, 1))
  rest <- substr(s, 2, nchar(s))
  return(paste0(first.letter, rest))
}

CapitalizeRemoveUnderscores <- function(s){
  # split string by underscore, capitalize first letter of each split word, 
  # paste together without underscore
  words <- strsplit(s, split = "_")[[1]]
  words.cap <- sapply(words, function(word) Capitalize(word))
  return(paste0(words.cap, collapse = " "))
}



# Load data ---------------------------------------------------------------

scale.factor <- 100
pseudocount <- 1

exprs.path <- "data/kallisto_encode/encode_rnaseq_abundances.annotated.sorted.merged.txt"
exprs.dat <- LoadEncodeRnaseq(exprs.path, long.format = TRUE, transform.log2 = FALSE)
exprs.signal <- LoadEncodeRnaseq(exprs.path, long.format = FALSE, transform.log2 = FALSE)
# exprs.dat <- read.table(exprs.path, header = TRUE)


# Sum across transcripts to get gene expression ---------------------------

exprs.bygene <- exprs.dat %>%
  group_by(gene, tissue) %>%
  summarise(tpm = sum(tpm))

exprs.bygene$logtpm <- log2(exprs.bygene$tpm * scale.factor + pseudocount)

# save Robj for easy retrieval later
# dat.encode <- exprs.bygene; save(dat.encode, file = "Robjs/dat.encode.Robj"); rm(dat.encode)


# Plot genes --------------------------------------------------------------

jgene <- "Onecut1"
jgene <- "Alb"
jgene <- "Ddc"
jgene <- "Hnf4a"
jgene <- "Mef2c"
jgene <- "Onecut2"
jgene <- "Nkx2-1"
jgene <- "Nr2f1"
jgene <- "Foxa2"
jgene <- "Rfxank"
# RFX1..5_RFXANK_RFXAP.p2

dat.sub <- SortByTissue(subset(exprs.bygene, gene == jgene & !tissue %in% c("Colon", "Duodenum", "Small Intestine", "Stomach", "Placenta")))

ggplot(dat.sub, aes(x = tissue, y = tpm)) + geom_bar(stat = "identity") + ggtitle(jgene) + xlab("Tissue") + ylab("TPM") +
  theme_gray(24) + theme(axis.text.x=element_text(angle=90,vjust = 0, hjust = 1))

# Quick PCA ---------------------------------------------------------------

s <- prcomp(exprs.signal, center = TRUE, scale. = FALSE)
plot(s$rotation[, 1], s$rotation[, 2])

screeplot(s)
textplot(s$rotation[, 1], s$rotation[, 2], rownames(s$rotation), cex=0.7)