# Load Kallisto abundance estimates and try to find first exons from there
# 2015-05-12

# detach("package:plyr", unload=TRUE)  # run this if you get dplyr bugs
library(dplyr)
library(mixtools)

# Functions ---------------------------------------------------------------

source("scripts/functions/GetTissueTimes.R")
source("scripts/functions/LoadArrayRnaSeq.R")
source("scripts/functions/FitRhythmic.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/AlternativeFirstExonsFunctions.R")
source("scripts/functions/MixtureModelFunctions.R")

LoadKallisto <- function(path.kallisto){
  source("scripts/functions/ConvertRNASeqTissueNamesToArray.R")
  if (missing(path.kallisto)){
    path.kallisto <- "data/alternative_exon_usage//abundance.merged.annotated.sorted.txt"
  }
  dat.kallisto <- read.table(path.kallisto, header = TRUE)
  
  # BEGIN: break matrix into a metadata and tpm estimates: makes converting to long easier
  dat.meta <- dat.kallisto[, c("chromo", "start", "end", "gene_name", "strand")]
  rownames(dat.meta) <- dat.kallisto$target_id
  
  colnames.meta <- c("chromo", "start", "end", "gene_name", "strand", "target_id")
  colnames.tpm <- colnames(dat.kallisto)[which(!colnames(dat.kallisto) %in% colnames.meta)]
  
  dat.tpm <- dat.kallisto[, colnames.tpm]
  rownames(dat.tpm) <- dat.kallisto$target_id
  # END: break matrix into a metadata and tpm estimates
  
  tissues <- sapply(colnames(dat.tpm), function(s) strsplit(s, '_')[[1]][[1]])
  tissues <- ConvertRNASeqTissueNamesToArray(tissues)
  times <- GetTimes(colnames(dat.tpm), get_unique=FALSE)
  
  tpm.long <- data.frame(transcript_id = rep(dat.kallisto$target_id, ncol(dat.tpm)),
                         chromo = rep(dat.kallisto$chromo, ncol(dat.tpm)),
                         start = rep(dat.kallisto$start, ncol(dat.tpm)), 
                         end = rep(dat.kallisto$end, ncol(dat.tpm)), 
                         gene_name = rep(dat.kallisto$gene_name, ncol(dat.tpm)),
                         strand = rep(dat.kallisto$strand, ncol(dat.tpm)),
                         tissue = rep(tissues, each = nrow(dat.tpm)),
                         time = as.numeric(rep(times, each = nrow(dat.tpm))),
                         tpm = unlist(dat.tpm))
  return(tpm.long)
}

IsTissueSpecific2 <- function(jdf, pval.min = 1e-5, pval.max = 0.05, cutoff = 4.89){
  # given list of pvals, check if it contains pvals less than pval.min and greater than pval.max.
  # check if pval contains values less than pval.min (rhythmic) and greater than pval.max (flat)
  # if so, return TRUE, otherwise FALSE
  avg.exprs <- jdf$int.rnaseq
  pvals <- jdf$pval
  pvals.filt <- pvals[which(!is.nan(pvals) & avg.exprs > cutoff)]
  if (min(pvals.filt) < pval.min & max(pvals.filt) > pval.max){
    jdf$is.tissue.spec.circ <- rep(TRUE, length(pvals))
  } else {
    jdf$is.tissue.spec.circ<- rep(FALSE, length(pvals))
  }
  return(jdf)
}

IsRhythmic2 <- function(pval, avg.exprs, pval.min = 1e-5, cutoff = 6){
  # ask if gene is rhythmic, not rhythmic or NA (lowly expressed or 0)
  # check if pval is less than pval.min
  if (is.nan(pval) | avg.exprs < cutoff){
    return(NA)
  }
  if (pval < pval.min){
    return(TRUE)
  } else {
    return(FALSE)
  }
}

ModelRhythmicity <- function(dat, jformula=tpm_normalized ~ is.rhythmic){
  # check it contains both Rhythmic and NotRhythmic elements
  if (length(unique(dat$is.rhythmic)) != 2){
    return(data.frame(int = NA, coef = NA, pval = NA))
  }
  jfit <- lm(data = dat, formula = jformula)
  # int <- summary(jfit)$mat["(Intercept)", "Coef"]
  # jcoef <- summary(jfit)$mat["rhythmic.or.notRhythmic", "Coef"]
  int <- coef(jfit)[[1]]
  jcoef <- coef(jfit)[[2]]
  pval <- summary(jfit)$coefficients["is.rhythmicTRUE", "Pr(>|t|)"]
  return(data.frame(int = int, coef = jcoef, pval = pval))
}

CalculateFractionIsoformUsage <- function(tpm, pseudocount = 0){
  # given tpm of a sample across known isoforms, compute
  # fraction of isoform usage. 
  # 
  # Add pseudo count to increase robustness to lowly expressed genes
  tpm_normalized <- (tpm + pseudocount) / (sum(tpm + pseudocount))
}

PlotDiagnostics <- function(jgene, jtranscript){
  # Plot diagnostic plots, given a gene and transcript
  # jgene <- "Sorbs1"; jtranscript="ENSMUST00000165469"  # Liver and BFAT looks rhythmic, but assigned as "not rhythmic"
  test.gene <- subset(tpm.afe.filt, gene_name == jgene)
  test <- subset(test.gene, transcript_id == jtranscript)
  
  print(PlotGeneAcrossTissues(subset(dat.long, gene == jgene)))
  
  print(PlotTpmAcrossTissues(test.gene, jtitle = paste(jgene, jtranscript), log2.transform=FALSE))
  print(PlotTpmAcrossTissues(test, jtitle = paste(jgene, jtranscript), log2.transform=TRUE))
  
  print(ggplot(test.gene, aes(y = tpm_normalized, x = transcript_id)) + 
    geom_boxplot() + 
    facet_wrap(~tissue) + 
    ggtitle(jgene))
  
  print(ggplot(test, 
         aes(y = tpm_normalized, x = is.rhythmic)) + 
    geom_boxplot() + 
    ggtitle(paste(jgene, jtranscript, "\n", unique(test$tissue[which(test$is.rhythmic == TRUE)]))))
}

# Load data ---------------------------------------------------------------

tpm.long <- LoadKallisto()
dat.long <- LoadArrayRnaSeq()

head(tpm.long)

jgene <- "Ddc"
test <- subset(tpm.long, gene_name == jgene)

test.afe <- test %>%
  group_by(time, tissue) %>%
  mutate(tpm_normalized = (tpm + 1) / sum(tpm + 1))
  
ggplot(test.afe, aes(y = tpm_normalized, x = transcript_id)) + geom_boxplot() + facet_wrap(~tissue) + ggtitle(jgene)


# Fit rhythmic ------------------------------------------------------------

dat.rhyth <- FitRhythmicDatLong(dat.long)

# dat.rhyth$is.rhythmic <- sapply(dat.rhyth$pval, IsRhythmic2)
dat.rhyth$is.rhythmic <- mapply(IsRhythmic2, pval = dat.rhyth$pval, avg.exprs = dat.rhyth$int.rnaseq)


# Find cutoff for rhythmicity ---------------------------------------------

# plot(density(dat.rhyth$int.rnaseq))  # visual inspection
# # mixture of two gaussians?
# mixmdl <- normalmixEM(dat.rhyth$int.rnaseq, lambda = c(0.50, 0.50), mu = c(2, 10), k = 2)
# plot(mixmdl,which=2)
# lines(density(dat.rhyth$int.rnaseq), lty=2, lwd=2)
# 
# cutoff <- optimize(ShannonEntropyMixMdl, interval = c(2, 8), mixmdl = mixmdl, maximum = TRUE)
# cutoff <- cutoff$maximum  # cutoff = 4.883356
# cutoff <- 4.883356
# 
# print(paste("Cutoff:", cutoff))

# Which genes are rhythmic in a tissue-specific manner? -------------------

dat.rhyth <- dat.rhyth %>%
  group_by(gene) %>%
  do(IsTissueSpecific2(., pval.min = 1e-5, pval.max = 0.05, cutoff = 6))

tissue.spec.circ.genes <- unique(dat.rhyth[which(dat.rhyth$is.tissue.spec.circ == TRUE), ]$gene)

print(paste("Number of tissue specific genes:", length(tissue.spec.circ.genes)))

# Calculate fractional isoform usage --------------------------------------


tpm.afe <- tpm.long %>%
  group_by(gene_name, tissue, time) %>%
  mutate(tpm_normalized = CalculateFractionIsoformUsage(tpm, pseudocount = 1))

jgene <- "Insig2"
ggplot(subset(tpm.afe, gene_name == jgene), aes(y = tpm_normalized, x = transcript_id)) + 
  geom_boxplot() + 
  facet_wrap(~tissue) + 
  ggtitle(jgene)


# Find associations between fractional isoform usage and rhythmici --------

# FIRST: assign genes as "rhythmic" or "not rhythmic"
tpm.afe.filt <- subset(tpm.afe, gene_name %in% tissue.spec.circ.genes)

rhythmic.dic <- setNames(object = dat.rhyth$is.rhythmic, nm = paste(dat.rhyth$gene, dat.rhyth$tissue, sep = ';'))

tpm.afe.filt$is.rhythmic <- rhythmic.dic[paste(tpm.afe.filt$gene_name, tpm.afe.filt$tissue, sep = ';')]
# DONE FIRST

# 2: Associate tissue-specific rhythms with fractional isoform usage
fit.afe <- tpm.afe.filt %>%
  filter(!is.na(is.rhythmic)) %>%  # filter out NA's. Onl consider TRUE and FALSE
  group_by(gene_name, transcript_id) %>%
  do(ModelRhythmicity(., jformula = tpm_normalized ~ is.rhythmic))

# 3. summarize by choosing the top for each gene
fit.afe.summary <- fit.afe %>%
  group_by(gene_name) %>%
  do(SubsetMinPval(jdf = .))

fit.afe.summary$pval.adj <- p.adjust(fit.afe.summary$pval)

head(data.frame(fit.afe.summary[order(fit.afe.summary$pval), ]), n = 100)

# How many make it past threshold?
pval.adj.thres <- 0.05

genes.tested <- unique(fit.afe.summary$gene_name)
n.hits <- length(which(fit.afe.summary$pval.adj <= pval.adj.thres))

sprintf("%s hits found out of %s tissue-specific circadian genes tested. %f percent", 
        n.hits, length(genes.tested), 100 * (n.hits / length(genes.tested)))


# Sanity checks -----------------------------------------------------------

# print top 100 hits
top.hits <- data.frame(fit.afe.summary[order(fit.afe.summary$pval), ])
top.hits <- head(top.hits, n = 10)

start <- Sys.time()
pdf("plots/alternative_exon_usage/kallisto_diagnostics4.pdf", width = 21, height = 7, paper = "USr")
jgenes <- as.character(top.hits$gene_name); jtranscripts <- as.character(top.hits$transcript_id)
mapply(PlotDiagnostics, jgene = jgenes, jtranscript = jtranscripts)
dev.off()
print(Sys.time() - start)

jgene <- "Abi2"
jgene <- "Bcat1"; jtranscript="ENSMUST00000123930"
jgene <- "Slc6a19"; jtranscript="ENSMUST00000124406"
jgene <- "Csrp3"; jtranscript="ENSMUST00000032658"
jgene <- "Insig2"; jtranscript="ENSMUST00000003818"
# jgene <- "Hnf4a"; jtranscript="ENSMUST00000109411"
jgene <- "Clcn1"; jtranscript="ENSMUST00000031894"
jgene <- "Ift80"; jtranscript="ENSMUST00000136448"
jgene <- "Il18r1"; jtranscript="ENSMUST00000087983"  # false positive
jgene <- "Atp6v1c2"; jtranscript="ENSMUST00000020884"
jgene <- "Rnf24"; jtranscript="ENSMUST00000110194"
jgene <- "Slc12a6"; jtranscript="ENSMUST00000053666"  # false positive? not really rhythmic this transcript
jgene <- "Insig2"; jtranscript="ENSMUST00000071064"  # probably real
jgene <- "Gne"; jtranscript="ENSMUST00000173274"  # probably real
jgene <- "Sparc"; jtranscript="ENSMUST00000125787"  # false positive
jgene <- "Sorbs1"; jtranscript="ENSMUST00000165469"  # Liver and BFAT looks rhythmic, but assigned as "not rhythmic"

# PlotTpmAcrossTissues(subset(tpm.afe, gene_name == jgene & transcript_id == jtranscript), log2.transform=FALSE)

# check rhythmicity status
subset(dat.rhyth, gene == jgene)

# jgene2 <- "Dbp"
# jgene2 <- "Elovl3"
# jgene2 <- "Rgs16"
# PlotTpmAcrossTissues(subset(tpm.afe, gene_name == jgene2), jtitle = jgene2, log2.transform=TRUE)
# PlotGeneAcrossTissues(subset(dat.long, gene == jgene2))

PlotGeneAcrossTissues(subset(dat.long, gene == jgene))

test.gene <- subset(tpm.afe.filt, gene_name == jgene)
test <- subset(test.gene, transcript_id == jtranscript)

PlotTpmAcrossTissues(test, jtitle = paste(jgene, jtranscript), log2.transform=TRUE)
PlotTpmAcrossTissues(test.gene, jtitle = paste(jgene, jtranscript), log2.transform=TRUE)
PlotTpmAcrossTissues(test.gene, jtitle = paste(jgene, jtranscript), log2.transform=FALSE)

ggplot(test.gene, aes(y = tpm_normalized, x = transcript_id)) + 
  geom_boxplot() + 
  facet_wrap(~tissue) + 
  ggtitle(jgene)

ggplot(test, aes(y = tpm_normalized, x = transcript_id)) + 
  geom_boxplot() + 
  facet_wrap(~tissue) + 
  ggtitle(paste(jgene, jtranscript))

ggplot(test, 
       aes(y = tpm_normalized, x = is.rhythmic)) + 
  geom_boxplot() + 
  ggtitle(paste(jgene, jtranscript, "\n", unique(test$tissue[which(test$is.rhythmic == TRUE)])))

ggplot(test.gene, aes(x = log2(tpm + 0.001))) + 
  geom_density() + 
  facet_wrap(~tissue) +
  ggtitle(paste(jgene))

# Plot density of TPM reads -----------------------------------------------

plot(density(log2(tpm.long$tpm + 0.01)))



