# 2017-05-01
# Jake Yeung
# Find target genes of MARA output

library(hash)
library(dplyr)
library(ggplot2)
library(plotrix)

source("scripts/functions/ListFunctions.R")
source("scripts/functions/HardcodedConstants.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/PlotFunctions.R")
source("/home/yeung/projects/sleep_deprivation/scripts/functions/DatabaseFunctions.R")

do.center <- TRUE

distfilt <- 40000
jweight <- 0
use.sql <- TRUE
jmod <- "Liver_SV129"
jcutoff <- 2
jcutoff.low <- 0

jcutoffstr <- paste(jcutoff, jcutoff.low, sep = ".")
jmodstr <- gsub(",", "-", jmod)

suffix <- GetSuffix(jweight, use.sql, jmodstr, jcutoffstr)
E.subdir <- GetESubDir(do.center, jmodstr, jweight)

maindir <- "/home/yeung/projects/tissue-specificity/Robjs/dhs_peaks"
prefix <- "liver.spec.peaks"

inf1 <- file.path(maindir, paste0(prefix, suffix, ".Robj"))

# Load matrces ------------------------------------------------------------

load("Robjs/liver_kidney_atger_nestle/dat.long.liverkidneyWTKO.bugfixed.Robj", v=T); dat.wtko <- dat.long; rm(dat.long)
load("Robjs/liver_kidney_atger_nestle/fits.bytiss.bugfixed.Robj", v=T)
fits.bytiss <- subset(fits.bytiss, tissue == "Liver_SV129" & gene != "")

# Load database -----------------------------------------------------------

# get tisspeaks and tissgenes
load(inf1, verbose = T)
liver.peaks <- unique(as.character(peaksgenes$peak))
liver.genes <- unique(as.character(peaksgenes$gene))

# distribution of phases of genes WITH liver peaks
ggplot(subset(fits.bytiss, gene %in% liver.genes), aes(x = phase)) + geom_histogram(bins = 30)

jcex <- 2
circular_phase24H_histogram(subset(fits.bytiss, gene %in% liver.genes)$phase, jtitle = paste(""), cex.axis = jcex, cex.lab = jcex, color_hist = "red")


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
print(paste("Collected", length(unique(as.character(N.long.filt$gene))), "genes and ", length(unique(as.character(N.long.filt$peak))), "peaks"))
print(str(N.long.filt))
print(Sys.time() - start)


# Merge and find top hits to explain MARA results -------------------------

N.sum <- N.long.filt %>%
  group_by(gene, motif) %>%
  summarise(sitecount = sum(sitecount))


# Associate target genes with phases --------------------------------------

phases <- hash(as.character(fits.bytiss$gene), fits.bytiss$phase)
means <- hash(as.character(fits.bytiss$gene), fits.bytiss$int)

N.sum$phase <- sapply(N.sum$gene, function(g) phases[[g]])
N.sum$mean <- sapply(N.sum$gene, function(g) means[[g]])

top.n <- 50


# Plot top hits -----------------------------------------------------------




# targets of ROR
jmotif <- "RORA.p2"
Nsub <- subset(N.sum, motif == jmotif) %>% 
  arrange(desc(sitecount))
ggplot(Nsub, aes(x = log2(sitecount))) + geom_histogram(bins = 20) + ggtitle(paste("Sitecounts of", jmotif))
target.genes <- Nsub$gene[1:top.n]
ggplot(subset(Nsub, gene %in% target.genes), aes(x = mean)) + geom_histogram(bins = 25) + ggtitle(paste("Means Motif:", jmotif))
circular_phase24H_histogram(subset(Nsub, gene %in% target.genes)$phase, jtitle = paste("Phase Motif", jmotif), cex.axis = jcex, cex.lab = jcex, color_hist = "red")

# targets of Ebox
jmotif <- "bHLH_family.p2"
Nsub <- subset(N.sum, motif == jmotif) %>% 
  arrange(desc(sitecount))
ggplot(Nsub, aes(x = log2(sitecount))) + geom_histogram(bins = 20) + ggtitle(paste("Sitecounts of", jmotif))
target.genes <- Nsub$gene[1:top.n]
ggplot(subset(Nsub, gene %in% target.genes), aes(x = mean)) + geom_histogram(bins = 25) + ggtitle(paste("Means Motif:", jmotif))
circular_phase24H_histogram(subset(Nsub, gene %in% target.genes)$phase, jtitle = paste("Phase Motif", jmotif), cex.axis = jcex, cex.lab = jcex, color_hist = "red")

# targets of AR
jmotif <- "AR.p2"
Nsub <- subset(N.sum, motif == jmotif) %>% 
  arrange(desc(sitecount))
ggplot(Nsub, aes(x = log2(sitecount))) + geom_histogram(bins = 20) + ggtitle(paste("Sitecounts of", jmotif))
target.genes <- Nsub$gene[1:top.n]
ggplot(subset(Nsub, gene %in% target.genes), aes(x = mean)) + geom_histogram(bins = 25) + ggtitle(paste("Means Motif:", jmotif))
circular_phase24H_histogram(subset(Nsub, gene %in% target.genes)$phase, jtitle = paste("Phase Motif", jmotif), cex.axis = jcex, cex.lab = jcex, color_hist = "red")

# targets of FOXA2
jmotif <- "FOXA2.p3"
Nsub <- subset(N.sum, motif == jmotif) %>% 
  arrange(desc(sitecount))
ggplot(Nsub, aes(x = log2(sitecount))) + geom_histogram(bins = 20) + ggtitle(paste("Sitecounts of", jmotif))
target.genes <- Nsub$gene[1:top.n]
ggplot(subset(Nsub, gene %in% target.genes), aes(x = mean)) + geom_histogram(bins = 25) + ggtitle(paste("Means Motif:", jmotif))
circular_phase24H_histogram(subset(Nsub, gene %in% target.genes)$phase, jtitle = paste("Phase Motif", jmotif), cex.axis = jcex, cex.lab = jcex, color_hist = "red")
