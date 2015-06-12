# Load Kallisto abundance estimates and try to find first exons from there
# 2015-05-12

# detach("package:plyr", unload=TRUE)  # run this if you get dplyr bugs
library(dplyr)
library(mixtools)

# Functions ---------------------------------------------------------------

source("scripts/functions/LoadKallisto.R")
source("scripts/functions/GetTissueTimes.R")
source("scripts/functions/LoadArrayRnaSeq.R")
source("scripts/functions/FitRhythmic.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/AlternativeFirstExonsFunctions.R")
source("scripts/functions/MixtureModelFunctions.R")

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
jgene <- "Srsf11"
# ggplot(subset(tpm.afe, gene_name == jgene), aes(y = tpm_normalized, x = transcript_id)) + 
#   geom_boxplot() + 
#   facet_wrap(~tissue) + 
#   ggtitle(jgene)

ggplot(subset(tpm.afe, gene_name == jgene), aes(y = tpm_normalized, x = tissue)) + 
  geom_boxplot() + 
  facet_wrap(~transcript_id) + 
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
top.hits <- top.hits[1:707, ]

start <- Sys.time()
pdf("plots/alternative_exon_usage/kallisto_diagnostics.top707.pdf", width = 21, height = 7, paper = "USr")
jgenes <- as.character(top.hits$gene_name); jtranscripts <- as.character(top.hits$transcript_id)
dat.tpm <- subset(tpm.afe.filt, gene_name %in% jgenes); dat.arrayrnaseq <- subset(dat.long, gene %in% jgenes)
mapply(PlotDiagnostics, jgene = jgenes, jtranscript = jtranscripts, 
       MoreArgs = list(dat.tpm = dat.tpm, dat.arrayrnaseq = dat.arrayrnaseq))
dev.off()
print(Sys.time() - start)


# Individual checks -------------------------------------------------------


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

# # PlotTpmAcrossTissues(subset(tpm.afe, gene_name == jgene & transcript_id == jtranscript), log2.transform=FALSE)
# 
# # check rhythmicity status
# subset(dat.rhyth, gene == jgene)
# 
# # jgene2 <- "Dbp"
# # jgene2 <- "Elovl3"
# # jgene2 <- "Rgs16"
# # PlotTpmAcrossTissues(subset(tpm.afe, gene_name == jgene2), jtitle = jgene2, log2.transform=TRUE)
# # PlotGeneAcrossTissues(subset(dat.long, gene == jgene2))
# 
# PlotGeneAcrossTissues(subset(dat.long, gene == jgene))
# 
# test.gene <- subset(tpm.afe.filt, gene_name == jgene)
# test <- subset(test.gene, transcript_id == jtranscript)
# 
# PlotTpmAcrossTissues(test, jtitle = paste(jgene, jtranscript), log2.transform=TRUE)
# PlotTpmAcrossTissues(test.gene, jtitle = paste(jgene, jtranscript), log2.transform=TRUE)
# PlotTpmAcrossTissues(test.gene, jtitle = paste(jgene, jtranscript), log2.transform=FALSE)
# 
# ggplot(test.gene, aes(y = tpm_normalized, x = transcript_id)) + 
#   geom_boxplot() + 
#   facet_wrap(~tissue) + 
#   ggtitle(jgene)
# 
# ggplot(test, aes(y = tpm_normalized, x = transcript_id)) + 
#   geom_boxplot() + 
#   facet_wrap(~tissue) + 
#   ggtitle(paste(jgene, jtranscript))
# 
# ggplot(test, 
#        aes(y = tpm_normalized, x = is.rhythmic)) + 
#   geom_boxplot() + 
#   ggtitle(paste(jgene, jtranscript, "\n", unique(test$tissue[which(test$is.rhythmic == TRUE)])))
# 
# ggplot(test.gene, aes(x = log2(tpm + 0.001))) + 
#   geom_density() + 
#   facet_wrap(~tissue) +
#   ggtitle(paste(jgene))
# 
# # Plot density of TPM reads -----------------------------------------------
# 
# plot(density(log2(tpm.long$tpm + 0.01)))



