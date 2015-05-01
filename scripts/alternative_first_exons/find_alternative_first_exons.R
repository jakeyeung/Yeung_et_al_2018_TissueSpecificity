# 2015-04-28
# find_alternative_first_exons.R

setwd("~/projects/tissue-specificity/")

library(ggplot2)
library(dplyr)
library(mixtools)
library(gplots)
library(parallel)
library(biglm)
# Functions ---------------------------------------------------------------

source("scripts/functions/LoadArrayRnaSeq.R")
source("scripts/functions/LoadAndHandleData.R")
source("scripts/functions/MixtureModelFunctions.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/MakeCluster.R")
source("scripts/functions/ReadListToVector.R")
source("scripts/functions/GrepRikGenes.R")
source("scripts/functions/FitRhythmic.R")


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
if (exists("dat.fitrhyth")) rm(dat.fitrhyth.split); rm(dat.long.by_genetiss); rm(dat.long.by_genetiss.split)  # GC

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


# Filter rhythmic genes with multistarts ----------------------------------

# cov.long.filts <- subset(cov.long, gene %in% common.genes)

# Filter lowly expressed genes --------------------------------------------

# cov.long.filt <- subset(cov.long.filts, rnaseq_reads >= 2^cutoff)  # cutoff established in log2 scale


# Optionaly do not filter anything ----------------------------------------

cov.long.filt <- cov.long


# Normalize exon reads ----------------------------------------------------

# cov.long.filt$reads_norm <- cov.long.filt$reads / cov.long.filt$rnaseq_reads  # naive

# n_starts tells us to filter this out because there are no alterantive exons
cov.normreads <- cov.long %>%
  group_by(tissue, gene, time) %>%
  mutate(norm_reads = Normalize(reads), n_starts = length(reads)) %>%
  filter(n_starts > 1)  # has no alternative first exons


# Find cutoff for background expression -----------------------------------

# this doesn't work as well unless we take the average across time...
# # find cutoff
# normreads.vec <- log2(cov.normreads$reads + 1)
# # takes ~1 minute
# mixmdl.normreads <- normalmixEM(normreads.vec, lambda = c(0.5, 0.5), mu = c(0.1, 6), k = 2)
# plot(mixmdl.normreads, which = 2)
# lines(density(normreads.vec), lty = 2, lwd = 2)
# cutoff.normreads <- optimize(ShannonEntropyMixMdl, interval = c(1, 5), mixmdl = mixmdl.normreads, maximum = TRUE)

# cutoff.normreads <- 2^(cutoff.normreads$maximum)

cutoff.normreads <- 2.019573
print(paste("mean reads cutoff:", cutoff.normreads))  # 2.01


# Remove lowly expressed genes --------------------------------------------

cov.normreads.filt <- cov.normreads %>%
  group_by(tissue, gene) %>%
  mutate(mean_reads.gene = mean(reads)) %>%
  filter(mean_reads.gene > cutoff.normreads)


# Filter only for genes of interest ---------------------------------------

cov.normreads.filt.common.genes <- subset(cov.normreads.filt, gene %in% common.genes)

# Ask if gene is rhythmic in that tissue ----------------------------------

dat.fitrhyth.filt <- FitDfToMatrix(dat.fitrhyth, common.genes)

rhythmic.or.not.mat <- apply(dat.fitrhyth.filt, 1, function(x) RhythmicOrNot(pval = x[1], amp = x[2]))

start <- Sys.time()
rnames <- paste(cov.normreads.filt.common.genes$tissue, cov.normreads.filt.common.genes$gene, sep = "-")
rhythmic.or.not <- rhythmic.or.not.mat[rnames]
print(Sys.time() - start)
if (exists(rhythmic.or.not)) rm(rhythmic.or.not.mat)

cov.normreads.filt.rhyth <- cbind(cov.normreads.filt.common.genes, rhythmic.or.not)


# Model reads on rhythmic or not ------------------------------------------

# test <- subset(cov.normreads.filt.rhyth, gene == "Ddc" & transcript == "ENSMUST00000178704" & !(is.na(rhythmic.or.not)))
# test <- subset(cov.normreads.filt.rhyth, gene == "Ddc" & transcript == "ENSMUST00000134121" & !(is.na(rhythmic.or.not)))
# test <- subset(cov.normreads.filt.rhyth, gene == "Cdh1" & transcript == "ENSMUST00000000312" & !(is.na(rhythmic.or.not)))
# test <- subset(cov.normreads.filt.rhyth, gene == "Adra1b" & transcript == "ENSMUST00000067258" & !(is.na(rhythmic.or.not)))
# test <- subset(cov.normreads.filt.rhyth, gene == "Csrp3" & transcript == "ENSMUST00000167786" & !(is.na(rhythmic.or.not)))
# test <- subset(cov.normreads.filt.rhyth, gene == "Arhgef10l" & transcript == "ENSMUST00000154979" & !(is.na(rhythmic.or.not)))

# for sanity checking
jgene <- "Asl"; jtrans <- "ENSMUST00000160557"
test <- subset(cov.normreads.filt.rhyth, gene == jgene & transcript == jtrans & !(is.na(rhythmic.or.not)))

PlotGeneAcrossTissues(subset(dat.long, gene == jgene))
ggplot(test, aes(x = rhythmic.or.not, y = norm_reads)) + geom_point() + ggtitle(paste(jgene, jtrans))
# ggplot(test, aes(x = rhythmic.or.not, y = log2(norm_reads))) + geom_point() + ggtitle(paste(jgene, jtrans))
fit.test <- biglm(formula = log2(norm_reads) ~ rhythmic.or.not, data = test)
(summary(fit.test))

# run model on full dataset
fit.afe <- cov.normreads.filt.rhyth %>%
  filter(!(is.na(rhythmic.or.not))) %>%
  group_by(transcript, gene) %>%
  do(FitRhythNonRhyth(jdf = .)) %>%
  filter(!is.na(pval))

# show top hits
(head(data.frame(fit.afe[order(fit.afe$pval), ]), n = 50))
(head(data.frame(fit.afe[order(fit.afe$coef, decreasing = TRUE), ]), n = 50))

# summarize by choosing the top for each gene
fit.afe.summary <- fit.afe %>%
  group_by(gene) %>%
  do(SubsetMinPval(jdf = .))

# plot histogram of pvalues
plot(density(fit.afe$pval))  # NICE

head(data.frame(fit.afe.summary[order(fit.afe.summary$pval), ]), n = 50)

# save(cov.normreads.filt.rhyth, file = "results/alternative_exon_usage/cov.normreads.filt.rhyth.Robj")

# How many peaks correlate with rhythmic genes? ---------------------------

pval.adj <- p.adjust(fit.afe.summary$pval)
fit.afe.summary$pval.adj <- pval.adj
(head(data.frame(fit.afe.summary[order(fit.afe.summary$pval.adj), ]), n = 50))
n.hits <- nrow(subset(fit.afe.summary, pval.adj < 0.05))
print(paste0("Number of hits:", n.hits))
sprintf("%s/%s are hits", n.hits, length(common.genes))

# # Calculate maximum difference --------------------------------------------
# 
# cov.normreads.sub <- subset(cov.normreads.filt, gene %in% common.genes)
# cov.normreads.by_gene <- group_by(cov.normreads.sub, gene)
# 
# # find AFEs by "minimum correlation"
# cov.mincor <- do(.data = cov.normreads.by_gene, GetMinCor(df = .))  # super slow
# 
# # Calculate log2 fold change ----------------------------------------------
# 
# # find AFEs by log2 fold change between "rhythmic" and "not rhythmic" genes
# cov.normreads.by_genetiss <- group_by(cov.normreads, gene, tissue)
# 
# 
# # ask if rhythmic 20 seconds
# rhythmic.or.not <- apply(cov.normreads.by_gene, 1, GetRhythmicOrNot, fitdf = dat.fitrhyth.filt)
# 
# # Append to df
# cov.normreads.by_gene.rhyth <- cbind(cov.normreads.by_gene, rhythmic.or.not)
# 
# # Fin
# 
# 
# 
# # Show distributions ------------------------------------------------------
# 
# (head(as.data.frame(cov.mincor[order(cov.mincor$min.cor), ]), n = 50))
# plot(density(cov.mincor$min.cor), xlim=c(0, 1))
# 
# 
# # Do sanity checks on examples --------------------------------------------
# 
# jgene <- "Ddc"
# jgene <- "Sgk2"
# jgene <- "Gas7"
# jgene <- "Insig2"
# # plot exprs
# PlotGeneAcrossTissues(subset(dat.long, gene == jgene))
# # check normreads
# jdf <- subset(cov.normreads.filt, gene == jgene)
# print(data.frame(jdf))
# # check out the matrix for correlations
# m <- acast(jdf, transcript ~ tissue, value.var = "norm_reads")
# # check which two tissues were called as "most different"
# best.cor <- LoopCor(m, show.which = TRUE)
# # check which transcript accounts for largest difference
# (m.sub <- m[, best.cor$tissues])
# # return transcript with highest difference
# m.diffs <- apply(m.sub, 1, function(x) abs(log2(x[1] / x[2])))
# (m.diffs.max <- m.diffs[which(m.diffs == max(m.diffs))])
# transcript.max <- names(m.diffs.max)
# # get chromosome location of this transcript
# GetLocationFromAnnotation(bed, gene_name = jgene, transcript_id = transcript.max)
# # Heatmap the matrix of correlations
# my_palette <- colorRampPalette(c("white", "black"))(n = 299)
# heatmap.2(m, density.info = "density", trace = "none", margins = c(8, 14), main = jgene, col = my_palette, cexRow=1.1)
# 
# MaxDiff <- 
# 
# # Calculate ShannonEntropy ------------------------------------------------
# 
# cov.entropy <- summarise(cov.normreads, entropy = ShannonEntropy(norm_reads))
# 
# # Remove NaNs -------------------------------------------------------------
# 
# cov.long.filt.entropy <- cov.long.filt.entropy[which(!is.nan(cov.long.filt.entropy$entropy)), ]
# 
# # Sort by entropy ---------------------------------------------------------
# 
# print(head(cov.long.filt.entropy[order(cov.long.filt.entropy$entropy), ], n = 50))
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
