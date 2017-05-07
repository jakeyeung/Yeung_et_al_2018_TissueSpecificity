# 2017-05-05
# Jake Yeung
# find pairs: first get cutoffs for each TF, and then use cutoff to find "hits"
# plot "hits" on genome-browser.

rm(list=ls())

library(hash)
library(dplyr)
library(ggplot2)
library(plotrix)

setwd("/home/yeung/projects/tissue-specificity")

source("scripts/functions/ListFunctions.R")
source("scripts/functions/HardcodedConstants.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/PlotFunctions.R")
source("scripts/functions/HandleMotifNames.R")
source("scripts/functions/LoadActivitiesLong.R")
source("scripts/functions/SvdFunctions.R")
source("scripts/functions/CosSineFunctions.R")
source("scripts/functions/PlotActivitiesFunctions.R")
source("scripts/functions/LiverKidneyFunctions.R")
source("scripts/functions/GetTFs.R")
source("scripts/functions/ProteomicsFunctions.R")
source("scripts/functions/FisherTestSitecounts.R")
source("/home/yeung/projects/sleep_deprivation/scripts/functions/DatabaseFunctions.R")


do.center <- TRUE

distfilt <- 40000
jweight <- 0
use.sql <- TRUE
jmod <- "Liver_SV129,Liver_BmalKO"
jmod <- "Liver_SV129"
jcutoff <- 3
jcutoff.low <- 0
incl.promoters <- FALSE

jcutoffstr <- paste(jcutoff, jcutoff.low, sep = ".")
jmodstr <- gsub(",", "-", jmod)
jtissgeno <- strsplit(jmod, ",")[[1]][[1]]
jtiss.nogeno <- strsplit(strsplit(jmod, ",")[[1]], "_")[[1]][[1]]

suffix <- GetSuffix(jweight, use.sql, jmodstr, jcutoffstr, incl.promoters)
E.subdir <- GetESubDir(do.center, jmodstr, jweight)

# maindir <- "/home/yeung/projects/tissue-specificity/Robjs/dhs_peaks"
maindir <- "/home/yeung/data/tissue_specificity/tissuepeaksgenes"
prefix <- "liver.spec.peaks"

jmotif <- "AR"
jmotif <- "FOXO1.3.4"
jmotif <- "DBP"
jmotif <- "bHLH_family"

inf1 <- file.path(maindir, paste0(prefix, suffix, ".Robj"))

# get tisspeaks and tissgenes
load(inf1, verbose = T)
liver.peaks <- unique(as.character(peaksgenes$peak))
liver.genes <- unique(as.character(peaksgenes$gene))


# Load database -----------------------------------------------------------

inf <- "/home/shared/sql_dbs/closestbed_multiple_genes.genomewide.merged.motifindexed.sqlite3"
motevo.tbl <- LoadDatabase(inf)

print("Getting genes from database")
start <- Sys.time()
N.sub.lst <- expandingList()
for (jgene in liver.genes){
  N.long.filt.query <- filter(motevo.tbl, gene == jgene)  # peaks are not indexed, so dont take them
  N.sub.tmp <- collect(N.long.filt.query, n = Inf)
  N.sub.lst$add(N.sub.tmp)
}
N.long.filt <- N.sub.lst$as.list()
N.long.filt <- bind_rows(N.long.filt)

rm(N.sub.tmp, N.sub.lst)  # worth it? 

# filter peaks after querying database
N.long.filt <- subset(N.long.filt, peak %in% liver.peaks & dist <= distfilt)
N.long.filt$motif <- sapply(N.long.filt$motif, function(m) make.names(RemoveP2Name(m)))

print(paste("Collected", length(unique(as.character(N.long.filt$gene))), "genes and ", length(unique(as.character(N.long.filt$peak))), "peaks"))
print(str(N.long.filt))
print(Sys.time() - start)

# Get flat peaks
suffix.flat <- GetSuffix(jweight, use.sql, "flat", jcutoffstr, incl.promoters)
inf.flat <- file.path(maindir, paste0(prefix, suffix.flat, ".Robj"))
load(inf.flat, v=T)
flat.peaks <- unique(as.character(peaksgenes$peak))
flat.genes <- unique(as.character(peaksgenes$gene))

print("Getting flat genes from database")
start <- Sys.time()
N.sub.lst <- expandingList()
for (jgene in flat.genes){
  N.long.filt.query <- filter(motevo.tbl, gene == jgene)  # peaks are not indexed, so dont take them
  N.sub.tmp <- collect(N.long.filt.query, n = Inf)
  N.sub.lst$add(N.sub.tmp)
}
N.long.flat <- N.sub.lst$as.list()
N.long.flat <- bind_rows(N.long.flat)
N.long.flat <- subset(N.long.flat, peak %in% flat.peaks & dist <= distfilt)
N.long.flat$motif <- sapply(N.long.flat$motif, function(m) make.names(RemoveP2Name(m)))

# save(N.long.filt, N.long.flat, file="Robjs/N.long.flat.from.sql.Robj")
print(Sys.time() - start)

N.merged <- bind_rows(N.long.filt, N.long.flat)

N.merged <- N.merged %>%
  group_by(motif, gene, peak) %>%
  summarise(sitecount = sum(sitecount))

# assign peak/gene to liver or flat
livflat.hash <- hash(c(liver.genes, flat.genes), c(rep("liver", length(liver.genes)), rep("flat", length(flat.genes))))

N.merged$model <- sapply(N.merged$gene, function(g) livflat.hash[[g]])

N.mat <- dcast(N.merged, formula = "gene + peak + model ~ motif", value.var = "sitecount", fill = 0)
N.merged.filled <- melt(N.mat, id.vars = c("gene", "peak", "model"), variable.name = "motif", value.name = "sitecount")

cutoffs <- seq(from = 0, to = 1, by = 0.1)
N.fisher <- RunFisherOnPromoters(N.merged.filled, foreground.models = "liver", background.models = "flat", cutoffs = cutoffs, return.full.df = FALSE)

# Cutoff at ~0.25
ggplot(subset(N.fisher, motif == "RORA"), aes(x = cutoff, y = -log10(p.value))) + geom_point() + geom_line()

# test on single motif
jsub <- subset(N.merged.filled, motif == "RORA")
test <- FisherTestSitecounts(jsub, 0.1, show.table = TRUE)


# # Do enrichment analysis by hypergeometric test ---------------------------
# 
# N.gene <- N.long.filt %>%
#   group_by(gene, motif) %>%
#   summarise(sitecount = sum(sitecount))
#   
#   
# cutoffs <- seq(from = 0.4, to = 0.8, by = 0.1)
# N.sub <- RunFisherOnPromoters(N.long, foreground.models = models.tw, background.models = models.flat, cutoffs = cutoffs)
# 
# 
# # Load matrces ------------------------------------------------------------
# 
# load("Robjs/liver_kidney_atger_nestle/dat.long.liverkidneyWTKO.bugfixed.Robj", v=T); dat.wtko <- dat.long; rm(dat.long)
# dat.wtko <- StaggeredTimepointsLivKid(dat.wtko)
# 
# load("Robjs/liver_kidney_atger_nestle/fits.bytiss.bugfixed.Robj", v=T)
# fits.bytiss <- subset(fits.bytiss, tissue == "Liver_SV129" & gene != "")
# 
# prot.long <- LoadProteomicsData()
# prot.long <- subset(prot.long, geno == "WT" & tissue == "Liver")



