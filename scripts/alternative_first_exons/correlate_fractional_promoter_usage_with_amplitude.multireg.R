# Jake Yeung
# 2015-10-29
# Normalize amplitude by Nr1d1 DO NOT adjust for noisy genes

setwd("/home/yeung/projects/tissue-specificity")
library(dplyr)
library(hash)
library(ggplot2)
library(reshape2)
library(biglm)

source("scripts/functions/LoadArrayRnaSeq.R")
source("scripts/functions/FitRhythmic.R")
source("scripts/functions/LoadKallisto.R")
source("scripts/functions/SvdFunctions.R")
source("scripts/functions/AlternativeFirstExonsFunctions.R")
source("scripts/functions/PlotGeneAcrossTissues.R")

filt.tiss <- c("WFAT")


# Function ----------------------------------------------------------------



# Load dat ----------------------------------------------------------------

load("Robjs/dat.long.fixed_rik_genes.Robj", verbose=T)
load("Robjs/dat.fit.Robj", verbose=T)
load("Robjs/alt_promoter_usage.mean_peak_amp.no_wfat.Robj", verbose = T)   # tpm.afe, tpm.avg.filt, tpm.fit


# Avg across time ---------------------------------------------------------


tpm.afe.avg <- tpm.afe %>%
  group_by(gene_name, transcript_id, tissue) %>%
  summarise(tpm_norm.avg = mean(tpm_normalized),
            tpm_norm.var = var(tpm_normalized),
            tpm.avg = mean(tpm),
            tpm.var = var(tpm))
tpm.afe.avg <- tpm.afe.avg[order(tpm.afe.avg$tpm_norm.var, decreasing = T), ]


# Canonical correlation ---------------------------------------------------

key <- paste(dat.fit$tissue, dat.fit$gene, sep = ",")
val <- dat.fit$amp
amp.dic <- hash(key, val)
mean.dic <- hash(key, dat.fit$int.rnaseq)
tpm.afe.avg$amp <- mapply(function(jgene, jtiss) amp.dic[[paste(jtiss, jgene, sep = ",")]],
                          as.character(tpm.afe.avg$gene_name), as.character(tpm.afe.avg$tissue))
tpm.afe.avg$mean <- mapply(function(jgene, jtiss) mean.dic[[paste(jtiss, jgene, sep = ",")]],
                          as.character(tpm.afe.avg$gene_name), as.character(tpm.afe.avg$tissue))

# count promoters
tpm.counts <- subset(tpm.afe.avg, tissue == "Adr") %>%
  group_by(tissue, gene_name) %>% 
  summarise(counts = length(transcript_id))
tpm.counts <- tpm.counts[order(tpm.counts$counts, decreasing = T), ]

# Genome-wide now ---------------------------------------------------------


counts.dic <- hash(as.character(tpm.counts$gene_name), tpm.counts$counts)
tpm.afe.avg$nprom <- sapply(as.character(tpm.afe.avg$gene_name), function(jgene) counts.dic[[jgene]])

start <- Sys.time()
tpm.mr <- subset(tpm.afe.avg, nprom > 1) %>%
  group_by(gene_name) %>%
  do(CorrelateAmpPromMulti(.))
save(tpm.mr, file = "Robjs/tpm.mr.weighted.Robj")
print(Sys.time() - start)
