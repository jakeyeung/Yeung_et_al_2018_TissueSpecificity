# 2015-04-28
# find_alternative_first_exons.R

setwd("~/projects/tissue-specificity/")

library(ggplot2)
library(dplyr)
library(mixtools)
library(gplots)
library(parallel)
library(biglm)
library(reshape2)
# Functions ---------------------------------------------------------------

source("scripts/functions/LoadArrayRnaSeq.R")
source("scripts/functions/LoadAndHandleData.R")
source("scripts/functions/MixtureModelFunctions.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/MakeCluster.R")
source("scripts/functions/ReadListToVector.R")
source("scripts/functions/GrepRikGenes.R")
source("scripts/functions/FitRhythmic.R")
source("scripts/functions/TissueSpecificRhythmicsFunctions.R")
source("scripts/functions/AlternativeFirstExonsFunctions.R")

# Load RNASeq ---------------------------------------------------

dat.long <- LoadArrayRnaSeq()
dat <- LoadRnaSeq()

save(dat.long, file = "")


# Find rhythmic genes -----------------------------------------------------

dat.rhyth <- FitRhythmicDatLong(dat.long)

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

# save(cov.long, file="results/alternative_exon_usage/cov.long.Robj")

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

# What is distribution of start sites?
cov.normreads.startsites <- cov.long %>%
  group_by(gene) %>%
  summarise(n_starts = length(transcript)) %>%
  mutate(n_starts_norm_samp = n_starts / 96) %>%  # 96 samples, so divide that
  select_("gene", "n_starts_norm_samp")

head(cov.normreads.startsites[order(cov.normreads.startsites$n_starts, decreasing = TRUE), ])

# filter for common.genes
canonicalgenes <- unique(dat.long$gene)
cov.normreads.startsites.common <- subset(cov.normreads.startsites, gene %in% canonicalgenes)

hist(cov.normreads.startsites.common$n_starts_norm_samp, breaks = 44, freq = TRUE, 
     main = paste0("Number of transcripts per gene (N=", 
                   length(cov.normreads.startsites.common$gene), ")"), 
     ylim = c(0, 3000))


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

test <- subset(cov.normreads.filt.rhyth, gene == "Ddc" & transcript == "ENSMUST00000178704" & !(is.na(rhythmic.or.not)))
# test <- subset(cov.normreads.filt.rhyth, gene == "Ddc" & transcript == "ENSMUST00000134121" & !(is.na(rhythmic.or.not)))
# test <- subset(cov.normreads.filt.rhyth, gene == "Cdh1" & transcript == "ENSMUST00000000312" & !(is.na(rhythmic.or.not)))
# test <- subset(cov.normreads.filt.rhyth, gene == "Adra1b" & transcript == "ENSMUST00000067258" & !(is.na(rhythmic.or.not)))
# test <- subset(cov.normreads.filt.rhyth, gene == "Csrp3" & transcript == "ENSMUST00000167786" & !(is.na(rhythmic.or.not)))
# test <- subset(cov.normreads.filt.rhyth, gene == "Arhgef10l" & transcript == "ENSMUST00000154979" & !(is.na(rhythmic.or.not)))

# for sanity checking
# jgene <- "Asl"; jtrans <- "ENSMUST00000160557"
jgene <- "Fbxl21"; jtrans <- "ENSMUST00000045428"
jgene <- "Eya1"; jtrans <- "ENSMUST00000027066"
test <- subset(cov.normreads.filt.rhyth, gene == jgene & transcript == jtrans & !(is.na(rhythmic.or.not)))

PlotGeneAcrossTissues(subset(dat.long, gene == jgene))
ggplot(test, aes(x = rhythmic.or.not, y = norm_reads)) + geom_point() + ggtitle(paste(jgene, jtrans))
# ggplot(test, aes(x = rhythmic.or.not, y = log2(norm_reads))) + geom_point() + ggtitle(paste(jgene, jtrans))
fit.test <- lm(formula = log2(norm_reads) ~ rhythmic.or.not, data = test)
(summary(fit.test))

# run model on full dataset
fit.afe <- cov.normreads.filt.rhyth %>%
  filter(!(is.na(rhythmic.or.not))) %>%
  group_by(transcript, gene) %>%
  do(FitRhythNonRhyth(jdf = .)) %>%
  filter(!is.na(pval))

cov.normreads.filt.rhyth.tmp <- filter(cov.normreads.filt.rhyth, !(is.na(rhythmic.or.not)))
cov.normreads.filt.rhyth.tmp <- group_by(cov.normreads.filt.rhyth.tmp, transcript, gene)
cov.normreads.filt.rhyth.tmp <- do(.data = cov.normreads.filt.rhyth.tmp, FitRhythNonRhyth(jdf = .))
cov.normreads.filt.rhyth.tmp <- filter(!is.na(pval))

# show top hits
(head(data.frame(fit.afe[order(fit.afe$pval), ]), n = 50))
(head(data.frame(fit.afe[order(fit.afe$coef, decreasing = TRUE), ]), n = 50))

# summarize by choosing the top for each gene
fit.afe.summary <- fit.afe %>%
  group_by(gene) %>%
  do(SubsetMinPval(jdf = .))

# plot histogram of pvalues
plot(density(fit.afe$pval))  # NICE

fit.afe.summary <- data.frame(fit.afe.summary[order(fit.afe.summary$pval), ])
head(fit.afe.summary, n = 50)

# save(cov.normreads.filt.rhyth, file = "results/alternative_exon_usage/cov.normreads.filt.rhyth.Robj")

# How many peaks correlate with rhythmic genes? ---------------------------

pval.adj.cutoff <- 0.05
pval.adj <- p.adjust(fit.afe.summary$pval)
fit.afe.summary$pval.adj <- pval.adj
(head(data.frame(fit.afe.summary[order(fit.afe.summary$pval.adj), ]), n = 50))
n.hits <- nrow(subset(fit.afe.summary, pval.adj < pval.adj.cutoff))
print(paste0("Number of hits:", n.hits))
sprintf("%s/%s are hits", n.hits, length(common.genes))


# Plot heatmap of examples ------------------------------------------------

jgene <- "Ddc"
cov.sub <- subset(cov.normreads.filt.rhyth, gene == jgene)
cov.sub.sum <- cov.sub %>%
  group_by(transcript, tissue) %>%
  summarise(mean_reads = mean(norm_reads))
cov.m <- acast(cov.sub.sum, transcript ~ tissue, value.var = "mean_reads")
my_palette <- colorRampPalette(c("white", "black"))(n = 299)
heatmap.2(cov.m, density.info = "histogram", trace = "none", margins = c(5, 14), main = jgene, col = my_palette, cexRow=1)


# Explore sitecount differences in hit -------------------------------------------------------------

source("scripts/functions/PlotSitecounts.R")
source("scripts/functions/LoadSitecounts.R")

N.list <- LoadSitecounts(gene_ids=FALSE)  # list of N and N.promoter
N.annot <- LoadEnsemblToPromoter()
# save(N.annot, file="Robjs/N.annot.Robj")
# rownames(N.annot) <- N.annot$ensemblid
N <- as.matrix(N.list$N)
N.promoter <- N.list$N.promoter
if (exists(x = "N") & exists(x = "N.promoter")) rm(N.list)
head(N)

N.long <- LoadSitecountsPromotersLong()
# save(N.long, file="Robjs/N.long.Robj")

jgene <- "Adra1b"
jgene <- "Csrp3"
jgene <- "Slc25a25"
jgene <- "Fpgs"
jgene <- "Insig2"
jgene <- "Ddc"
par(mar=c(5.1, 12, 4.1, 2.1))

# print(N[grep(jgene, rownames(N)), ])
# print(N.promoter[grep(jgene, rownames(N)), ])

m <- N[grep(jgene, rownames(N)), ]
m <- m[, which(colSums(m) > 0.25)]
my_palette <- colorRampPalette(c("white", "black"))(n = 299)
heatmap.2(t(m), density.info = "density", trace = "none", margins = c(1, 14), main = jgene, col = my_palette, cexRow=0.75)


# Correlate alternative transcript usage to promoter sitecounts -----------

# Filter only for genes that contain a correlation between promoter usage
# and rhythmicity

genes.hit <- subset(fit.afe.summary, pval.adj < pval.adj.cutoff)$gene
max.pval <- max(subset(fit.afe.summary, pval.adj < pval.adj.cutoff)$pval)  # any pval below this is a hit
fit.afe.hits <- subset(fit.afe, gene %in% genes.hit)

fit.afe.hits$hit.or.not <- factor(mapply(HitOrNot, fit.afe.hits$coef, fit.afe.hits$pval, 
                                         MoreArgs = list(max.pval=max.pval, min.pval=max.pval)), 
                                  levels = c("NotHit", "Neg", "Pos"))

fit.afe.hits.filt <- subset(fit.afe.hits, hit.or.not %in% c("Neg", "Pos"))

# Annotate N.long to hit or not (subset makes things easier)
N.annot.sub <- subset(N.annot, ensemblid %in% fit.afe.hits.filt$transcript)
annot.dic <- setNames(fit.afe.hits.filt$hit.or.not, fit.afe.hits.filt$transcript)
N.annot.sub$hit.or.not <- sapply(N.annot.sub$ensemblid, function(x){
  return(annot.dic[[x]])
})
annot.dic <- setNames(N.annot.sub$hit.or.not, N.annot.sub$saeedid)
N.long.sub <- subset(N.long, promoterid %in% N.annot.sub$saeedid)
N.long.sub$hit.or.not <- sapply(N.long.sub$promoterid, function(x){
  x <- as.character(x)
  return(annot.dic[[x]])
})

# Fit linear model of motif and hit or not

# BEGIN TESTING MOTIFS INDIVIDUALLY
jmotif <- "RORA.p2"
jmotif <- "TFDP1.p2"
jmotif <- "NFIL3.p2"
jmotif <- "bHLH_family.p2"
jmotif <- "MTF1.p2"
jmotif <- "PAX8.p2"
jmotif <- "POU1F1.p2"
jmotif <- "HNF1A.p2"
test <- subset(N.long.sub, motif == jmotif)

ks.test <- KsTestPosNeg(jdf = test)
fit.test <- FitPosNeg(jdf = test)
ks.test; fit.test

ggplot(data = test, aes(x = hit.or.not, y = sitecount)) + 
  geom_boxplot() + 
  ggtitle(jmotif)
ggplot(data = test, aes(x = sitecount, colour = hit.or.not, fill = hit.or.not)) + 
  geom_density(alpha = 0.5) +
  ggtitle(jmotif)
# END TESTING MOTIFS INDIVIDUALLY

ks.motif <- N.long.sub %>%
  group_by(motif) %>%
  do(KsTestPosNeg(jdf = .)) %>%
  .[order(.$pval), ]
ks.motif$pval.adj <- p.adjust(ks.motif$pval)
head(ks.motif)

fit.motif <- N.long.sub %>%
  group_by(motif) %>%
  do(FitPosNeg(jdf = .)) %>%
  .[order(.$pval), ]
fit.motif$pval.adj <- p.adjust(fit.motif$pval)
head(fit.motif)

