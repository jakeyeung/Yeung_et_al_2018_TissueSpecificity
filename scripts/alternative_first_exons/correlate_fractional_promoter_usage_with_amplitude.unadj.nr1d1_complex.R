# Jake Yeung
# 2015-10-29
# Normalize amplitude by Nr1d1 DO NOT adjust for noisy genes

setwd("/home/yeung/projects/tissue-specificity")
library(dplyr)
library(hash)
library(ggplot2)
library(reshape2)

source("scripts/functions/LoadArrayRnaSeq.R")
source("scripts/functions/FitRhythmic.R")
source("scripts/functions/LoadKallisto.R")
source("scripts/functions/SvdFunctions.R")
source("scripts/functions/AlternativeFirstExonsFunctions.R")
source("scripts/functions/PlotGeneAcrossTissues.R")

filt.tiss <- c("WFAT")

# Get data ----------------------------------------------------------------


dat.long <- LoadArrayRnaSeq()
dat.fit <- FitRhythmicDatLong(dat.long)
tpm.merged <- LoadKallisto(path.kallisto = "data/alternative_exon_usage/abundance.merged.annotated.sorted.pythonmerged.bed")

# Get complex matrix ------------------------------------------------------
ref.gene <- "Nr1d1"
dat.fit.relamp <- GetRelampByGene(dat.fit, by.gene = ref.gene)
dat.fit.relamp <- dat.fit.relamp %>%
  group_by(gene) %>%
  mutate(max.exprs = max(int.rnaseq), min.pval = min(pval), max.relamp = max(relamp), relamp.norm = relamp / sum(relamp), amp.norm = amp / sum(amp))

genes.exprs <- unique(subset(dat.fit.relamp, max.exprs >= 4)$gene)
dat.complex <- TemporalToFrequencyDatLong(subset(dat.long, gene %in% genes.exprs), period = 24, n = 8, interval = 6, add.entropy.method = "array")

# normalize to reference gene (Nr1d1)
ref.amps <- subset(dat.complex, gene == ref.gene, select = c(tissue, gene, exprs.transformed))
ref.amps.dic <- hash(ref.amps$tissue, ref.amps$exprs.transformed)

dat.complex$exprs.transformed.norm <- mapply(function(tiss, exprs) exprs / ref.amps.dic[[tiss]], as.character(dat.complex$tissue), dat.complex$exprs.transformed)
dat.complex$mod.exprs.transformed <- Mod(dat.complex$exprs.transformed)
dat.complex$mod.exprs.transformed.norm <- Mod(dat.complex$exprs.transformed.norm)


# Annotate complex mat with mean value ------------------------------------

means.dic.keys <- paste(dat.fit.relamp$tissue, dat.fit.relamp$gene, sep = ";")
means.dic.vals <- dat.fit.relamp$int.rnaseq
means.dic <- hash(means.dic.keys, means.dic.vals)

dat.complex$mean <- mapply(function(gene, tissue) means.dic[[paste(c(tissue, gene), collapse = ";")]], 
                           as.character(dat.complex$gene), as.character(dat.complex$tissue))

# Get only genes where at least one tissue is rhythmic --------------------

rhyth.genes <- unique(subset(dat.complex, mod.exprs.transformed.norm >= 0.1)$gene)


# Calculate frac usage ----------------------------------------------------


tpm.afe <- tpm.merged %>%
  filter(gene_name %in% rhyth.genes) %>%
  group_by(gene_name, tissue, time) %>%
  mutate(tpm_normalized = CalculateFractionIsoformUsage(tpm, pseudocount = 1))

tpm.avg <- tpm.afe %>%
  group_by(gene_name, transcript_id, tissue) %>%
  summarise(tpm_norm.avg = mean(tpm_normalized), tpm_norm.var = var(tpm_normalized)) %>%
  filter(tpm_norm.avg < 1 & tpm_norm.var > 0)


# Annotate rel amp and mean exprs -----------------------------------------


# annotate with normalized relative amplitude and mean expression
relamp.dic.keys <- paste(dat.complex$tissue, dat.complex$gene, sep = ";")
relamp.dic.vals <- 2 * dat.complex$mod.exprs.transformed
relamp.dic <- hash(relamp.dic.keys, relamp.dic.vals)

tpm.avg$relamp <- mapply(function(gene, tissue) relamp.dic[[paste(c(tissue, gene), collapse = ";")]], 
                              as.character(tpm.avg$gene_name), as.character(tpm.avg$tissue))
tpm.avg$int.rnaseq <- mapply(function(gene, tissue) means.dic[[paste(c(tissue, gene), collapse = ";")]],
                       as.character(tpm.avg$gene_name), as.character(tpm.avg$tissue))

# Model frac usage to amp -------------------------------------------------

tpm.avg.filt <- subset(tpm.avg, int.rnaseq > 3 & !tissue %in% filt.tiss)
tpm.fit <- tpm.avg.filt %>%
  group_by(gene_name, transcript_id) %>%
  do(FitPromoterUsageToAmplitude(.))


# Show hits ---------------------------------------------------------------

tpm.summary <- tpm.fit %>%
  group_by(gene_name) %>%
  do(SubsetMinPval(jdf = .))

tpm.summary$pval.adj <- p.adjust(tpm.summary$pval)

pval.cutoff <- 0.005
tpm.summary.filt <- subset(tpm.summary, pval <= pval.cutoff)
sig.hits <- tpm.summary.filt$gene_name

head(data.frame(tpm.summary[order(tpm.summary$pval), which(colnames(tpm.summary) != "transcript_id")]), n = 100)
head(data.frame(tpm.summary[order(tpm.summary$relamp.range, decreasing = TRUE), which(colnames(tpm.summary) != "transcript_id")]), n = 100)
head(data.frame(tpm.summary[order(tpm.summary$tpm_norm.range, decreasing = TRUE), which(colnames(tpm.summary) != "transcript_id")]), n = 100)
head(data.frame(tpm.summary[order(tpm.summary$tpm_norm.range, decreasing = TRUE), ]), n = 100)


# Hits --------------------------------------------------------------------

plot(density(tpm.summary$pval[which(!is.na(tpm.summary$pval))]))

hits <- head(data.frame(tpm.summary[order(tpm.summary$pval), which(colnames(tpm.summary) != "transcript_id")]), n = 200)$gene_name

# pdf("plots/alternative_exon_usage/linear_model_complex_mat_top_hits.no_WFAT.real_amp.pdf")
# tpm.hits <- subset(tpm.afe, gene_name %in% hits)
# tpm.avg.filt.hits <- subset(tpm.avg.filt, gene_name %in% hits)
# dat.long.hits <- subset(dat.long, gene %in% hits)
# for (jgene in hits){
#   tpm.sub <- subset(tpm.fit.hits, gene_name == jgene)
#   jtranscript <- tpm.sub[order(abs(tpm.sub$relamp), decreasing = TRUE), ]$transcript_id[1]
#   PlotDiagnostics2(subset(tpm.hits, gene_name == jgene), 
#                    subset(tpm.avg.filt.hits, transcript_id == jtranscript), 
#                    subset(dat.long.hits, gene == jgene), jgene, jtranscript)
# }
# dev.off()

jgene <- "Upp2"
jgene <- "Ddc"
jgene <- "Insig2"
jgene <- "Tpm3"
jgene <- "Prkd3"
jgene <- "Cstad"
jgene <- "Npas2"

hits.sub <- hits[1:20]
tpm.hits <- subset(tpm.fit, gene_name %in% hits.sub)
dat.hits <- subset(dat.long, gene %in% hits.sub)
tpm.afe.hits <- subset(tpm.afe, gene_name %in% hits.sub)
tpm.avg.filt.hits <- subset(tpm.avg.filt, gene_name %in% hits.sub)
pdf("plots/alternative_exon_usage/linear_model_complex_mat_top_hits.no_WFAT.real_amp.pdf")
for (jgene in hits.sub){
  print(jgene)
  tpm.sub <- subset(tpm.fit, gene_name == jgene)
  jtranscript <- tpm.sub[order(abs(tpm.sub$relamp), decreasing = TRUE), ]$transcript_id[1]
  PlotDiagnostics2(tpm.afe, tpm.avg.filt, dat.long, jgene, jtranscript)  
}
dev.off()

jgene <- "Bpifa1"
PlotGeneAcrossTissues(subset(dat.long, gene == jgene))

# Plot hits in UCSC -------------------------------------------------------

# jtranscripts <- vector(mode = "numeric", length = length(hits.sub))
# for (i in 1:length(hits.sub)){
#   jgene <- hits.sub[i]
#   tpm.sub <- subset(tpm.fit, gene_name == jgene)
#   jtranscripts[i] <- as.character(tpm.sub[order(abs(tpm.sub$relamp), decreasing = TRUE), ]$transcript_id[1])
#   PlotGeneAcrossTissues(subset(dat.long, gene == jgene))
#   print(PlotTpmAcrossTissues(subset(tpm.afe, gene_name == jgene), jtitle = paste("log-scale", jgene, jtranscript, sep = "\n"), log2.transform=TRUE))
# }
# 
# plot_bed <- subset(tpm.merged, transcript_id %in% jtranscripts & tissue == "Liver" & time == 22, select = c("chromo", "start", "end", "gene_name"))
# bedToUCSC(plot_bed, "plots/alternative_exon_usage/hits_ucsc.pdf", leftwindow = 20000, rightwindow = 19999)


# Visualize TPM usage as orthogonal axes ----------------------------------

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

tpm.afe.avg$amp <- mapply(function(jgene, jtiss) amp.dic[[paste(jtiss, jgene, sep = ",")]],
                          as.character(tpm.afe.avg$gene_name), as.character(tpm.afe.avg$tissue))

# count promoters
tpm.counts <- subset(tpm.afe.avg, tissue == "Adr") %>%
  group_by(tissue, gene_name) %>% 
  summarise(counts = length(transcript_id))
tpm.counts <- tpm.counts[order(tpm.counts$counts, decreasing = T), ]

DoCanCor <- function(tpm.sub, xvar = "tpm_norm.avg", yvar = "amp"){
  # do canonical correlation of two matrices to find interesting
  # alternative promoter usage
  if (nrow(tpm.sub) == 0){
    return(data.frame(Cor = NA))
  }
  X <- dcast(data = tpm.sub, formula = transcript_id ~ tissue, value.var = xvar)
  rownames(X) <- X$transcript_id; X$transcript_id <- NULL
  jtrans <- tpm.sub$transcript_id[[1]]
  Y <- data.frame(subset(tpm.sub, transcript_id == jtrans, select = c("amp", "tissue")))
  rownames(Y) <- as.character(Y$tissue); Y$tissue <- NULL
  jcor <- cancor(tpm.afe.mat, dat.fit.mat, xcenter = F, ycenter = F)
  return(data.frame(Cor = jcor$cor))
}

tpm.afe.avg.split <- split(tpm.afe.avg, tpm.afe.avg$gene_name)

library(parallel)
tpm.cancor <- mclapply(tpm.afe.avg.split, function(x) DoCanCor(x), mc.cores = 30)
head(tpm.cancor)
save(tpm.cancor, file = "Robjs/tpm.cancor.Robj")

# # canonical correlation
# jgene <- "Insig2"
# for (jgene in head(tpm.counts)$gene_name){
#   # set up AFE usage matrix
#   tpm.afe.mat <- subset(tpm.afe.avg, gene_name == jgene)
#   tpm.afe.mat <- dcast(data = tpm.afe.mat, formula = transcript_id ~ tissue, value.var = "tpm_norm.avg")
#   rownames(tpm.afe.mat) <- tpm.afe.mat$transcript_id; tpm.afe.mat$transcript_id <- NULL
#   tpm.afe.mat <- t(tpm.afe.mat)  # n by p 
#   
#   dim(tpm.afe.mat)
#   #   if more promoters than samples, then do SVD to reduce the dimensions into a square.
#   if (ncol(tpm.afe.mat) >= 3){
#     print("Reducing dimensions...")
#     s <- svd(tpm.afe.mat)
#     rownames(s$u) <- rownames(tpm.afe.mat)
#   #     tpm.afe.mat <- s$u[, 1:nrow(tpm.afe.mat) - 1]
#       tpm.afe.mat <- s$u[, 1:3]
#   }
#   
#   # set up amplitude and mean
#   dat.fit.mat <- data.frame(subset(dat.fit, gene == jgene, select = c("amp", "int.rnaseq", "tissue")))
#   rownames(dat.fit.mat) <- as.character(dat.fit.mat$tissue); dat.fit.mat$tissue <- NULL
#   dat.fit.mat <- subset(dat.fit.mat, select = c("amp"))
#   #   print(dim(tpm.afe.mat))
#   #   print(dim(dat.fit.mat))
#   tpm.afe.mat <- as.matrix(tpm.afe.mat)
#   dat.fit.mat <- as.matrix(dat.fit.mat)
#   #   tpm.afe.mat <- sweep(tpm.afe.mat, MARGIN = 2, STATS = colMeans(tpm.afe.mat))
#   jcor <- cancor(tpm.afe.mat, dat.fit.mat, xcenter = F, ycenter = F)
#   
#   # check correlation
#   xcor <- tpm.afe.mat %*% jcor$xcoef
#   ycor <- dat.fit.mat %*% jcor$ycoef
#   # ycor <- as.matrix(dat.fit.mat)
#   #   cor(xcor, ycor)
#   #   dim(tpm.afe.mat)
#   print(cor(xcor, ycor))
#   print(jcor$cor)
#   print(PlotGeneAcrossTissues(subset(dat.long, gene == jgene)))
#   plot(xcor, ycor, main = paste(jgene, jcor$cor))
#   text(xcor, ycor, labels = rownames(xcor)) 
# }


