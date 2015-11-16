# 2015-11-10
# cluster_by_promoter_space.R
# Cluster tissues by promoter space

library(hash)
library(dplyr)
library(ggplot2)
library(reshape2)
library(mvtnorm)

setwd("/home/yeung/projects/tissue-specificity")

tiss <- c('Adr', 'WFAT', 'Lung','BS','Kidney','BFAT','Mus','Cere','Hypo','Aorta','Liver','Heart')  
# keep WFAT to get the dic without <NULL>. Remove WFAT using subsets later.


# Functions ---------------------------------------------------------------

source("scripts/functions/GetClockGenes.R")
source('scripts/functions/AlternativeFirstExonsFunctions.R')
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/NcondsFunctions.R")

GetGeneModelKeys <- function(dat, tiss){
  rhyth.tiss <- gsub(pattern = ";", replacement = ",", x = dat$model)
  rhyth.tiss <- strsplit(rhyth.tiss, ",")[[1]]
  
  is.rhyth <- tiss %in% rhyth.tiss
  return(data.frame(tissue = tiss, is.rhyth = is.rhyth))
}

# Load --------------------------------------------------------------------

load("Robjs/dat.long.fixed_rik_genes.Robj", verbose=T)
load("Robjs/tpm.afe.avg.Robj", verbose=T)
load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.15.phase_sd_maxdiff_avg.Robj", verbose=T)


# Needs a tougher amp filer -----------------------------------------------

fits.best$model <- mapply(FilterModelByAmp, fits.best$model, fits.best$param.list, MoreArgs = list(amp.cutoff = 0.2))
fits.best$n.params <- sapply(fits.best$model, function(m) return(length(strsplit(as.character(m), ";")[[1]])))
fits.best$n.rhyth <- sapply(fits.best$model, GetNrhythFromModel)
fits.best$amp.avg <- mapply(GetAvgAmpFromParams, fits.best$param.list, fits.best$model)
fits.best$phase.sd <- mapply(GetSdPhaseFromParams, fits.best$param.list, fits.best$model)
fits.best$phase.maxdiff <- mapply(GetMaxPhaseDiffFromParams, fits.best$param.list, fits.best$model)
fits.best$phase.avg <- mapply(GePhaseFromParams, fits.best$param.list, fits.best$model)

# Annotate tpm.afe.avg by model -------------------------------------------

# hash it
keys.df <- fits.best %>%
  group_by(gene) %>%
  do(GetGeneModelKeys(., tiss))
keys <- paste(keys.df$tissue, keys.df$gene, sep = ",")
vals <- as.numeric(keys.df$is.rhyth)
is.rhyth.dic <- hash(keys, vals)

# annotate it
tpm.afe.avg <- subset(tpm.afe.avg, gene_name %in% unique(fits.best$gene))
tpm.afe.avg$amp <- mapply(function(jtiss, jgene) is.rhyth.dic[[paste(jtiss, jgene, sep = ",")]], 
                          as.character(tpm.afe.avg$tissue), as.character(tpm.afe.avg$gene_name))



# Genome-wide k-means and Gaussian clustering -----------------------------

start <- Sys.time()
cutoff <- 3
tpm.fuzzy <- subset(tpm.afe.avg, nprom > 1 & tissue != "WFAT" & mean > cutoff) %>%
  group_by(gene_name) %>%
  do(RunFuzzyDistance(.))
save(tpm.fuzzy, file = "Robjs/tpm.fuzzy.bic_models.Robj")

tpm.gauss <- subset(tpm.afe.avg, nprom > 1 & tissue != "WFAT" & mean > cutoff) %>%
  group_by(gene_name) %>%
  do(sigs = CalculateGaussianCenters(.))
tpm.gauss2 <- subset(tpm.gauss, !is.na(sigs)) %>%
  group_by(gene_name) %>%
  do(CalculateGaussianDists(.))

genelist <- tpm.gauss2$gene_name
tpm.gauss <- cbind(tpm.gauss2, subset(tpm.gauss, gene_name %in% genelist, select = -gene_name))
save(tpm.gauss, file = "Robjs/tpm.gauss.bic_models.Robj")
print(Sys.time() - start)

# # Test out the k-means and Gaussian clustering ----------------------------
# 
# jgene <- "Ddc"
# jgene <- "Pkp4"
# jgene <- "Slc45a3"
# 
# PlotGeneAcrossTissues(subset(dat.long, gene == jgene))
# RunFuzzyDistance(subset(tpm.afe.avg, gene_name == jgene))
# out1 <- CalculateGaussianCenters(subset(tpm.afe.avg, gene_name == jgene))
# out2 <- CalculateGaussianDists(out1) 
# 
# cgenes <- GetClockGenes()
# cgenes <- c(cgenes, "Ddc", "Pkp4", "Slc45a3")
# fuzzy <- subset(tpm.afe.avg, gene_name %in% cgenes & nprom > 1 & tissue != "WFAT" & mean > 3.5) %>%
#   group_by(gene_name) %>%
#   do(RunFuzzyDistance(.))
# 
# gauss <- subset(tpm.afe.avg, gene_name %in% cgenes & nprom > 1 & tissue != "WFAT" & mean > 3.5) %>%
#   group_by(gene_name) %>%
#   do(sigs = CalculateGaussianCenters(.))
# gauss2 <- subset(gauss, !is.na(sigs)) %>%
#   group_by(gene_name) %>%
#   do(CalculateGaussianDists(.))
  
