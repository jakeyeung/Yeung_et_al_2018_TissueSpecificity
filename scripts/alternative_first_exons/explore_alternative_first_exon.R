# Afte rrunning find_alternative_first_exons.R, I saved my objects so we can explore them later.


# Functions ---------------------------------------------------------------

setwd("~/projects/tissue-specificity/")

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

library(ggplot2)
library(dplyr)
library(mixtools)
library(gplots)
library(parallel)
library(biglm)
library(reshape2)
# Functions ---------------------------------------------------------------



# Load my data ------------------------------------------------------------


load(file = "results/alternative_exon_usage/cov.normreads.filt.rhyth.Robj")  # cov.normreads.filt.rhyth
load(file = "docs/2015-03-05_lab_meeting/robjs/arrayrnaseq.Robj")  # mydat
dat.long <- mydat; rm(mydat)
load(file = "results/alternative_exon_usage//cov.long.Robj")
load(file = "results/alternative_exon_usage/N_and_N.promoter.Robj")


# Explore -----------------------------------------------------------------

head(cov.normreads.filt.rhyth)

# run model on full dataset
fit.afe <- cov.normreads.filt.rhyth %>%
  filter(!(is.na(rhythmic.or.not))) %>%
  group_by(transcript, gene) %>%
  do(FitRhythNonRhyth(jdf = .)) %>%
  filter(!is.na(pval))

# show top hits
(head(data.frame(fit.afe[order(fit.afe$pval), ]), n = 50))
# (head(data.frame(fit.afe[order(fit.afe$coef, decreasing = TRUE), ]), n = 50))

# summarize by choosing the top for each gene
fit.afe.summary <- fit.afe %>%
  group_by(gene) %>%
  do(SubsetMinPval(jdf = .))

# plot histogram of pvalues
plot(density(fit.afe$pval))  # NICE

