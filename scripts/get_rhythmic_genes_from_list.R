# plot_exprs_from_genelist.R
# March 3 2015
library(ggplot2)
library(plyr)
library(foreach)
library(doParallel)
setwd("/home/yeung/projects/tissue-specificity")


# Constants ---------------------------------------------------------------

T <- 24  # hours
omega <- 2 * pi / T
outfile <- "results/tissue_specific_rhythmic_genes/tissue_specific_rhythmic_genes.txt"

# Functions ---------------------------------------------------------------

source("scripts/functions/LoadArrayRnaSeq.R")
source("scripts/functions/FitRhythmic.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/SvdFunctions.R")
source("scripts/functions/ConvertLongToWide.R")
source("scripts/functions/PlotActivitiesFunctions.R")
source("scripts/functions/ReadListToVector.R")
source("scripts/functions/TissueSpecificRhythmicsFunctions.R")  # for Is.TissueSpecific or Is.RhythmicAcrossTissues
source("scripts/functions/MakeCluster.R")

# Set cluster -------------------------------------------------------------

MakeCluster()

# Load data ---------------------------------------------------------------


scripts.dir <- "scripts"
funcs.dir <- "functions"

genes <- ReadListToVector("data/gene_lists/genes_associated_with_multiple_promoters.txt")
# genes <- ReadListToVector("data/gene_lists/genes_with_multiple_promoters_aftermm10liftOver.txt")
dat <- LoadArrayRnaSeq()
dat.sub <- subset(dat, gene %in% genes)

# Fit Rhythmic ------------------------------------------------------------

# Fit rhythmic
sprintf("Running with %s cores.", ncores)
t.start <- Sys.time()
dat.sub.split <- split(dat.sub, dat.sub$tissue)
# remove WFAT
dat.sub.split$WFAT <- NULL
dat.sub.fit <- lapply(dat.sub.split, function(df){
  ddply(df, .(gene), FitRhythmic, .parallel=TRUE)
})
print(t.start - Sys.time())

# Create long dataframe
dat.sub.fit.long <- do.call(rbind, dat.sub.fit)

# Summarize each gene by average and SD of pval
cutoff.pval <- 1e-5
cutoff.amp <- 0.5

# Filter out neural tissues
# Filter out WFAT because it's a weird one.
# dat.sub.fit.long.filt <- subset(dat.sub.fit.long, ! tissue %in% c("BS", "Cere", "Hypo", "WFAT"))
dat.sub.fit.long.filt <- dat.sub.fit.long  # no filter

# save as object to load for presentations
save(dat.sub.fit.long, file = "results/FitRhythmicOutputs/FitRhythmic.dat.sub.fit.long.Robj")


# Filter for genes containing significant pvals AND non significant across tissues
tissue_specific_genes_or_not <- ddply(dat.sub.fit.long.filt, .(gene), summarise,
                                      IsTissueSpecific = Is.TissueSpecific(pval, amp))

tissue_specific_genes <- subset(tissue_specific_genes_or_not, IsTissueSpecific == TRUE)$gene
non_specific_genes <-  subset(tissue_specific_genes_or_not, IsTissueSpecific == FALSE)$gene

# write tissue_specific_genes to file
dat.sub.tiss <- subset(dat.sub.fit.long.filt, gene %in% tissue_specific_genes)
dat.sub.tiss <- dat.sub.tiss[order(dat.sub.tiss$pval), ]
write.table(x = dat.sub.tiss, file = outfile, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

rhythmic_genes <- ddply(dat.sub.fit.long.filt, .(gene), summarise,
                        IsRhythmic = Is.RhythmicAcrossTissues(pval, amp))


# Plot some genes and look for alternative promoter usage on the b --------

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



# # Doing SVD on this matrix. What happens? ---------------------------------
# 
# dat.sub.filt <- subset(dat.sub, gene %in% tissue_specific_genes)
# dat.sub.split <- split(dat.sub.filt, dat.sub.filt$tissue)
# 
# dat.split.proj <- lapply(dat.sub.split, function(x){ 
#   ddply(x, .(gene), ProjectToFrequency2, omega = omega, add.tissue = TRUE, .parallel = TRUE)
# })
# dat.proj <- do.call(rbind, dat.split.proj)
# 
# 
# 
# # Make wide format --------------------------------------------------------
# 
# dat.proj.wide <- ConvertLongToWide(dat.proj)
# 
# 
# 
# # Perform SVD -------------------------------------------------------------
# 
# tissues <- colnames(dat.proj.wide)
# genes <- rownames(dat.proj.wide)
# 
# s <- svd(dat.proj.wide)
# rownames(s$u) <- genes
# rownames(s$v) <- tissues
# 
# plot(s$d^2 / sum(s$d ^ 2), type = 'o')
# 
# for (comp in seq(11)){
#   PlotEigengene(s, comp)
#   PlotEigensamp(s, comp)
# }


