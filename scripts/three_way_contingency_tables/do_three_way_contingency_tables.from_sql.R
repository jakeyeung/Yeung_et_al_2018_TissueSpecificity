# 2016-04-21
# Jake Yeung
# pLDA was done before with cutoff at 0.5, but no real reason to do that. Redo with 0.1!
# Rerun from SQL database where FastaN is cleaned to prevent any NNNN or lower case in fasta file 

rm(list=ls())

# Ks <- c(200, 300)
# # Ks <- c(100, 150, 200, 250, 350, 400, 450, 500, 550, 600)
Ks <- c(300)
# Ks <- seq(100, 600, by = 50)
jmc.cores <- length(Ks)

args <- commandArgs(trailingOnly=TRUE)

# Parameters to tune
# nb.levels <- 3  # liver-rhyth, nonliver, liver-nonrhyth
# nb.levels <- 2  # liver-rhyth, liver-nonrhyth
# weight.cutoff <- 0.8  # makes cleaner??

nb.levels <- as.numeric(args[[1]])
weight.cutoff <- as.numeric(args[[2]])

# weight.cutoff <- 0  # more genes?
null.model <- "JI"  # JI or CI
flat.livamp.cutoff <- Inf  # gene must be VERY flat!

RunPoissonModel.JI <- function(dat){
  mod1 <- glm(freq ~ model + motif1 + motif2 + motif1 * motif2, data=dat, family=poisson())
  chisq.pval <- 1 - pchisq(mod1$deviance, mod1$df.residual)
  return(data.frame(deviance = mod1$deviance, pval = chisq.pval))
}

RunPoissonModel.CI <- function(dat){
  mod1 <- glm(freq ~ model + motif1 + motif2 + model * motif1 + model * motif2, data=dat, family=poisson())
  chisq.pval <- 1 - pchisq(mod1$deviance, mod1$df.residual)
  return(data.frame(deviance = mod1$deviance, pval = chisq.pval))
}

library(data.table)
setwd("/home/yeung/projects/tissue-specificity")

args <- commandArgs(trailingOnly=TRUE)
distfilt <- as.numeric(args[1])
jcutoff <- as.numeric(args[2])
jcutoff.low <- as.numeric(args[3])
distfilt <- 40000
jcutoff <- 3  # arbitrary
jcutoff.low <- 0  # arbitrary
# weight.cutoff <- 0
# jcutoff <- 2  # arbitrary
# jcutoff <- 3  # arbitrary
cleanup <- FALSE
writepeaks <- FALSE
jmethod <- "g=1001"
# jmethod <- "BIC"

merge.peaks <- FALSE


if (null.model == "JI"){
  RunPoissonModel <- RunPoissonModel.JI
} else if (null.model == "CI"){
  RunPoissonModel <- RunPoissonModel.CI
}

print(paste("Null model:", null.model))
print(RunPoissonModel)

if (is.na(distfilt)) stop("Distfilt must be numeric")

print(paste("Distance filter:", distfilt))
print(paste("DHS signal zscore cutoff:", jcutoff))

saveplot <- FALSE
saverobj <- FALSE
outdir <- "plots/penalized_lda/liver_kidney_wtko"
dir.create(outdir)
# outf <- paste0(outdir, "2D.posterior.multigene.distfilt.liverWTKO.", distfilt, ".cutoff.", jcutoff, ".cutofflow", jcutoff.low, ".method.", jmethod, ".pdf")
amp.min <- 0

# jmodels <- c("Kidney_SV129")
jmodels <- c("Liver_SV129")
if (jmodels == "Kidney_SV129"){
  rhyth.tiss <- c("Kidney")
  flat.tiss <- c("Liver")
} else if (jmodels == "Liver_SV129"){
  rhyth.tiss <- c("Liver")
  flat.tiss <- c("Kidney")
}
# if (saverobj){
#   outfile.robj <- paste0("Robjs/liver_kidney_atger_nestle/penalized_lda_mats.posterior.model.", jmodels[[1]], ".distfilt.", distfilt, ".", paste(rhyth.tiss, sep = "_"), ".cutoff", jcutoff, ".cutofflow", jcutoff.low, ".method.", jmethod, ".Robj")
#   if (file.exists(outfile.robj)) stop(paste0(outfile.robj, " exists. Exiting"))
# }
# 
library(MASS)
library(dplyr)
library(reshape2)
library(ggplot2)
library(hash)
library(parallel)
# library(penalizedLDA)

source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/PlotFunctions.R")
source("scripts/functions/ShannonEntropy.R")
source("scripts/functions/SitecountsFunctions.R")
source("scripts/functions/LdaFunctions.R")
source("scripts/functions/RemoveP2Name.R")
source("scripts/functions/PlotUCSC.R")  # CoordToBed
source("scripts/functions/GetTopMotifs.R")
source("scripts/functions/HandleMotifNames.R")
source("scripts/functions/ListFunctions.R")
source("scripts/functions/PlotActivitiesFunctions.R")
source("scripts/functions/LoadActivitiesLong.R")
source("scripts/functions/PlotUCSC.R")
source("scripts/functions/LiverKidneyFunctions.R")
source("/home/yeung/projects/sleep_deprivation/scripts/functions/DatabaseFunctions.R")


# Functions ---------------------------------------------------------------

colMax <- function(dat){
  return(apply(dat, 1, max))
}

# Load --------------------------------------------------------------------

# save_N_on_posterior_cutoff_0.1.R saves Robj image. Here we laod it up
# Do Penalized LDA as before  ---------------------------------------------

# from multigene_analysis.play_with_parameters.R 
if (!exists("fits.best")){
  load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.15.phase_sd_maxdiff_avg.Robj", v=T)
  fits.best.hog <- fits.best
  
  load("Robjs/liver_kidney_atger_nestle/fits.long.multimethod.filtbest.staggeredtimepts.bugfixed.Robj", v=T)
  fits.best.orig <- fits.long.filt
  fits.best <- fits.long.filt; rm(fits.long.filt)
  fits.best <- subset(fits.best, method == jmethod)
  fits.best <- subset(fits.best, weight >= weight.cutoff)
} 
if (!exists("dat.long")){
  load("Robjs/dat.long.fixed_rik_genes.Robj", v=T)
  dat.long.hog <- dat.long
  load("Robjs/liver_kidney_atger_nestle/dat.long.liverkidneyWTKO.bugfixed.Robj", v=T)
  dat.orig <- dat.long
  dat.long <- CollapseTissueGeno(dat.long, keep.tissue.col = TRUE)
  dat.long <- StaggeredTimepointsLivKid(dat.long)
} 
if (!exists("S.long")) load("Robjs/S.long.multigene.filt.50000.Robj", v=T)
if (!exists("N.long.filt")){
  liver.flat.genes <- unique(as.character(subset(fits.best, model %in% c("", "Liver_SV129"))$gene))
  inf <- "/home/shared/sql_dbs/closestbed_multiple_genes.genomewide.merged.motifindexed.sqlite3"
  motevo.tbl <- LoadDatabase(inf)
  
  print("Getting genes from database")
  start <- Sys.time()
  N.sub.lst <- expandingList()
  for (jgene in liver.flat.genes){
    N.long.filt.query <- filter(motevo.tbl, gene == jgene)  # peaks are not indexed, so dont take them
    N.sub.tmp <- collect(N.long.filt.query, n = Inf)
    N.sub.lst$add(N.sub.tmp)
  }
  N.long.filt <- N.sub.lst$as.list()
  N.long.filt <- bind_rows(N.long.filt)
  
}
# get fit from f24

load("Robjs/liver_kidney_atger_nestle/fits.bytiss.bugfixed.Robj", v=T)
fits.bytiss <- subset(fits.bytiss, tissue %in% "Liver_SV129" & gene != "")
liver.amp <- hash(as.character(fits.bytiss$gene), fits.bytiss$amp)
liver.phase <- hash(as.character(fits.bytiss$gene), fits.bytiss$phase)

print("Assigning livamp to fitsbet")
fits.best$livamp <- sapply(as.character(fits.best$gene), function(g){
  jamp <- liver.amp[[g]]
  if (is.null(jamp)){
    jamp <- NA
  }
  return(jamp)  
})
print(head(fits.best))

# Get genes and peaks -----------------------------------------------------


jgenes <- as.character(subset(fits.best, model %in% jmodels)$gene)
jgenes.flat <- as.character(subset(fits.best, model == "" & livamp <= flat.livamp.cutoff)$gene)


print(paste("Rhythmic genes:", length(jgenes)))
print(paste("Flat genes:", length(jgenes.flat)))


# print(paste("Rhythmic genes:", length(jgenes)))
# print(paste("Flat genes:", length(jgenes.flat)))

# get peaks near gene
S.sub <- subset(S.long, gene %in% jgenes & dist < distfilt)
S.sub.flat <- subset(S.long, gene %in% jgenes.flat & dist < distfilt)
jpeaks <- as.character(unique(S.sub)$peak)  # 192726 peaks for Liver genes within 50 kb away
jpeaks.flat <- as.character(unique(S.sub.flat$peak))

print(paste("number of peaks in liver-specific rhythmic genes", length(jpeaks)))
print(paste("number of peaks in flat genes", length(jpeaks.flat)))

N.sub <- subset(N.long.filt, gene %in% jgenes & peak %in% jpeaks)  # should not have any missing peaks
N.sub.flat <- subset(N.long.filt, gene %in% jgenes.flat & peak %in% jpeaks.flat)

# Clean up ram ------------------------------------------------------------
if (cleanup){
  rm(S.long, N.long.filt)
}



# Identify tissue-specific peaks ------------------------------------------


# take peaks with Liver signal greater than cutoff
jtiss <- levels(S.sub$tissue)
tiss.i <- which(jtiss %in% rhyth.tiss)
others.i <- which(jtiss %in% flat.tiss)

S.sub.liverpeaks <- S.sub %>%
  group_by(peak, gene) %>%  # tissue order as "Cere", "Heart", "Kidney", "Liver", "Lung", "Mus"
  filter(min(zscore[tiss.i]) > jcutoff & max(zscore[others.i]) < jcutoff.low)

S.sub.nonliverpeaks <- S.sub %>%
  group_by(peak, gene) %>%
  # filter(min(zscore[others.i]) > jcutoff & max(zscore[tiss.i]) > jcutoff)
  # filter(min(zscore[others.i]) > 0.5 * jcutoff & max(zscore[tiss.i]) < jcutoff.low)
  filter(min(zscore[others.i]) > 0.5 * jcutoff & max(zscore[tiss.i]) < jcutoff.low)
  # filter(min(zscore[others.i]) > jcutoff)

jtiss.flat <- levels(S.sub.flat$tissue)
if (identical(jtiss, jtiss.flat) == FALSE){
  print("This shouldnt be necessary")
  tiss.i <- which(jtiss %in% rhyth.tiss)
  others.i <- which(jtiss %in% rhyth.tiss)
}

S.sub.flat.liverpeaks <- S.sub.flat %>%
  group_by(peak, gene) %>%
  filter(min(zscore[tiss.i]) > jcutoff & max(zscore[others.i]) < jcutoff.low)
  # filter(min(zscore[tiss.i]) > jcutoff)

S.sub.flat.nonliverpeaks <- S.sub.flat %>%
  group_by(peak, gene) %>%
  filter(min(zscore[others.i]) > jcutoff & max(zscore[tiss.i]) < jcutoff.low)
  # filter(min(zscore[others.i]) > jcutoff)


# Merge peaks if necessary ------------------------------------------------



# Mreg should contain peaks -----------------------------------------------

jpeak <- "chr1:72173321-72173821"
subset(S.long, peak == jpeak & tissue %in% c("Liver", "Kidney"))
subset(S.sub, peak == jpeak & tissue %in% c("Liver", "Kidney"))


# Do 3-way contingency tables ---------------------------------------------




rhyth.peaks <- as.character(unique(S.sub.liverpeaks$peak))
flat.peaks <- as.character(unique(S.sub.flat.liverpeaks$peak))
kidney.peaks <- as.character(unique(S.sub.nonliverpeaks$peak))
peaks.all <- c(rhyth.peaks, flat.peaks, kidney.peaks)
# remove duplicate peaks
peaks.all <- peaks.all[!duplicated(peaks.all)]

# rhythmic peaks first
N.sub.liver <- subset(N.sub, peak %in% rhyth.peaks)
N.sub.nonliver <- subset(N.sub, peak %in% kidney.peaks)
N.sub.flat.sub <- subset(N.sub.flat, peak %in% flat.peaks)

N.sub.liver$model <- "rhyth"
N.sub.nonliver$model <- "kidney"
N.sub.flat.sub$model <- "flat"

if (nb.levels == 3){
  N.merged <- rbind(N.sub.liver, N.sub.nonliver, N.sub.flat.sub)
} else if (nb.levels == 2){
  N.merged <- rbind(N.sub.liver, N.sub.flat.sub)
}
# N.merged <- rbind(subset(N.sub, peak %in% peaks.all), subset(N.sub.flat, peak %in% peaks.all))

# Optionally merge peaks, and add motifs?
if (merge.peaks){
  N.merged <- N.merged %>%
    group_by(motif, gene, model) %>%
    summarise(sitecount = sum(sitecount)) %>%
    mutate(peak = gene)
}

N.merged <- N.merged %>%
  group_by(peak, model, motif, gene) %>%
  summarise(sitecount = sum(sitecount)) %>%
  group_by(motif) %>%
  arrange(desc(sitecount)) %>%
  mutate(motif.rank = seq(length(sitecount)))

# take liver peaks and flat peaks
# K <- 200  # take top 1000 top guys
mclapply(Ks, function(K){
  
  top.bot <- c("atop", "zbottom")
  N.merged$in.list <- sapply(N.merged$motif.rank, function(m) ifelse(m < K, "atop", "zbottom"))
  
  start <- Sys.time()
  
  source("scripts/functions/GetTFs.R")
  N.merged$motif <- sapply(as.character(N.merged$motif), function(m) RemoveP2AndCommaBracesDashes(m))
  
  motif.pairs <- combn(as.character(unique(N.merged$motif)), m = 2, simplify = FALSE)
  
  # calculate atop, zbottom 2 by 2 by N tables. Where N is number of models. 
  N.mat.all <- dcast(N.merged, formula = gene + peak + model ~ motif, value.var = "in.list", fill = top.bot[2])
  
  start <- Sys.time()
  N.mat.freqs <- lapply(motif.pairs, function(motif.pair){
    motif.pair.str <- paste(motif.pair, collapse = ";")
    N.mat.freq <- N.mat.all %>%
      group_by_(.dots = c("model", motif.pair)) %>%
      summarise(freq = length(gene)) %>%
      mutate(pair = motif.pair.str)
    # expect 8 rows: atop/zbottom, zbottom/atop, atop/atop/ zbottom/zbottom for rhyth and flat
    # expect 12 rows if nb.levels = 3  (4 per level)
    if (nrow(N.mat.freq) == 4 * nb.levels){
      return(N.mat.freq)
    } else {
      warning("N.mat.freq must be 8 or 12 rows!")
      return(data.frame())
    }
  })
  print(Sys.time() - start)
  
  N.mat.freqs <- rbindlist(N.mat.freqs)
  colnames(N.mat.freqs) <- c("model", "motif1", "motif2", "freq", "pair")
  
  fits <- N.mat.freqs %>%
    group_by(pair) %>%
    do(RunPoissonModel(.))
  
  save(fits, N.mat.all, N.mat.freqs, file = paste0("Robjs/three_way_cooccurence/three.way.cooccurrence.fastabugfixed.UseSql.nmodels.", nb.levels, ".K.", K, ".weight.", weight.cutoff, ".MergePeaks.", merge.peaks, ".nullmodel.", null.model, ".flatampmax.", flat.livamp.cutoff,  ".withNmatallNmatfreqs.RemoveZeroCounts.Robj"))
}, mc.cores = jmc.cores)

# tf <- "RORA"
# tf <- "DBP"
# tf <- "SRF"
# tf <- "FOXA2"
# tf <- "CUX2"
# tf <- "RORA"
# tf <- "bHLH_family"
# tf <- "FOXA2"
# tf <- sample(unique(N.merged$motif), 1)
# tf <- "ONECUT"
# x <- data.frame(head(fits[grepl(tf, fits$pair), ] %>% arrange(desc(deviance)), n = 100))
# x$padj <- p.adjust(x$pval)

# # What is the co-occurrence table for hit pairs?
# 
# 
# 
# 
# # motif.pair <- c("HNF4A.p2", "RORA.p2") 
# # motif.pair <- c("HNF4A_NR2F1,2.p2", "RORA.p2")
# # motif.pair <- c("ONECUT1,2.p2", "RORA.p2")
# fits <- lapply(motif.pairs, function(motif.pair){
#   N.mat <- dcast(subset(N.merged, motif %in% motif.pair), formula = gene + peak + model ~ motif, value.var = "in.list", fill = top.bot[2])
#   colnames(N.mat) <- c("gene", "peak", "model", "motif1", "motif2")
#   
#   # do linear model
#   N.mat.freq <- N.mat %>%
#     group_by(model, motif1, motif2) %>%
#     summarise(freq = length(gene))
#   mod1 <- glm(freq ~ model + motif1 + motif2 + motif1 * motif2, data=N.mat.freq, family=poisson())
#   chisq.pval <- 1 - pchisq(mod1$deviance, mod1$df.residual)
#   return(data.frame(pair = paste(motif.pair, collapse = ";"), 
#                     deviance = mod1$deviance,
#                     pval = chisq.pval))
# })
# print(Sys.time() - start)
# 
# dir.create("Robjs/three_way_cooccurence")
# save(fits, file = "Robjs/three_way_cooccurence/three.way.cooccurrence.again.Robj")
# 
# # 
# # 
# # # N.mat <- dcast(subset(N.merged, motif %in% motif.pair), formula = gene + peak + model ~ motif, value.var = "in.list")
# # 
# # N.table <- xtabs(~motif1 + motif2 + model, N.mat)  # full model
# # N.table <- xtabs(~motif1 + motif2, N.mat)  # null model
# # 
# # debug <- TRUE
# # if (debug){
# #   # this is 2by2 to test with heart data: http://data.princeton.edu/wws509/notes/c5.pdf
# #   # N.table <- N.table[, , 1]
# #   # dat.byhand <- data.frame(serum = rep(c("low", "high"), each = 2), 
# #                            # disease = rep(c("yes", "no")),
# #                            # freq = c(51, 992, 41, 245))
# #   # N.table <- xtabs(freq ~ serum + disease, dat.byhand)
# #   # N.table <- matrix(c(51, 992, 41, 245), nrow = 2, ncol = 2, byrow = TRUE)
# #   
# #   # # this is 3by3 test
# #   dat.byhand <- data.frame(social = rep(c("lower", "lower.mid", "upper.mid", "higher"), each = 4),
# #                            parental = rep(c("low", "high"), each = 2),
# #                            college.plans = rep(c("no", "yes")),
# #                            freq = c(749, 35, 233, 133, 627, 38, 330, 303, 420, 37, 374, 467, 153, 26, 266, 800))
# #   N.table <- xtabs(freq ~ social + parental + college.plans, dat.byhand)
# #   # N.table <- xtabs(freq ~ parental + college.plans + social, dat.byhand)
# # }
# # 
# # # test <- loglm(~motif1 + motif2 + model + motif1 * motif2, N.table)
# # # test <- loglm(~motif1 + motif2 + model, N.table)
# # 
# # # calculate by hand: mutual independence
# # obs.table <- N.table / sum(N.table)  # easy
# # # expected table for mutual independence, do it margin by margin (3 times)
# # # exp.table <- N.table
# # # for (i in 1:2){
# # #   for (j in 1:2){
# # #       exp.table[i, j] <- sum(obs.table[i,]) * sum(obs.table[, j])
# # #   }
# # # }
# # 
# # exp.table <- N.table
# # for (i in 1:dim(N.table)[1]){
# #   for (j in 1:dim(N.table)[2]){
# #     for (k in 1:dim(N.table)[3]){
# #       # exp.table[i, j, k] <- (sum(N.table[i, ,]) / sum(N.table)) * (sum(N.table[, j, ]) / sum(N.table)) * (sum(N.table[, , k]) / sum(N.table))
# #       exp.table[i, j, k] <- sum(obs.table[i, ,]) * sum(obs.table[, j, ]) * sum(obs.table[, , k])
# #     }
# #   }
# # }
# # freq.table <- sum(N.table) * exp.table
# # # compare to expected
# # # freq.table.check <- sum(N.table[, 1]) * sum(N.table[1, ]) / sum(N.table)
# # # freq.table.check2 <- N.table[1, 1] * exp.table[1, 1]
# # 
# # logL.obs <- sum(N.table * log(obs.table))
# # logL.exp <- sum(N.table * log(exp.table))
# # deviance <- 2 * (logL.obs - logL.exp)  # ok
# # logL.ratio <- 2 * sum(N.table * log(N.table / freq.table))
# # chisq.pval <- 1 - pchisq(deviance, prod(dim(N.table)) - (sum(dim(N.table))) + 2)
# # 
# # 
# # # Easier way: poisson models ----------------------------------------------
# # # http://stats.stackexchange.com/questions/8342/appropriate-way-to-deal-with-a-3-level-contingency-table
# # # http://data.princeton.edu/wws509/notes/c5.pdf
# # 
# # # get frequency 
# # N.mat.freq <- N.mat %>%
# #   group_by(model, motif1, motif2) %>%
# #   summarise(freq = length(gene))
# # 
# # # mod0 <- glm(freq ~ model + motif1 + motif2, data=N.mat.freq, family=poisson())
# # mod1 <- glm(freq ~ model + motif1 + motif2 + motif1 * motif2, data=N.mat.freq, family=poisson())
# # 
# # chisq.pval <- 1 - pchisq(mod1$deviance, mod1$df.residual)
# # 
# # N.table.motifs <- xtabs(~motif1 + motif2, N.mat)
# # N.table.model <- xtabs(~model, N.mat)
