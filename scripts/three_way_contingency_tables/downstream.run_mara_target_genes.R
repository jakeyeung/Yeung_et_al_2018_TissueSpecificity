# 2017-05-01
# Jake Yeung
# Find target genes of MARA output

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
source("/home/yeung/projects/sleep_deprivation/scripts/functions/DatabaseFunctions.R")


do.center <- TRUE

distfilt <- 40000
jweight <- 0
use.sql <- TRUE
jmod <- "Liver_SV129,Liver_BmalKO"
jmod <- "Liver_SV129"
jcutoff <- 1.5
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

# Load matrces ------------------------------------------------------------

load("Robjs/liver_kidney_atger_nestle/dat.long.liverkidneyWTKO.bugfixed.Robj", v=T); dat.wtko <- dat.long; rm(dat.long)
dat.wtko <- StaggeredTimepointsLivKid(dat.wtko)

load("Robjs/liver_kidney_atger_nestle/fits.bytiss.bugfixed.Robj", v=T)
fits.bytiss <- subset(fits.bytiss, tissue == "Liver_SV129" & gene != "")

prot.long <- LoadProteomicsData()
prot.long <- subset(prot.long, geno == "WT" & tissue == "Liver")

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
N.long.filt$motif <- sapply(N.long.filt$motif, function(m) make.names(RemoveP2Name(m)))

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


# Any mean exprs bias between kidney and liver? ---------------------------

dat.mean <- dat.wtko %>%
  group_by(tissue, geno, gene) %>%
  summarise(exprs = mean(exprs)) %>%
  mutate(tissgeno = paste(tissue, geno, sep = "_"))

ggplot(subset(dat.mean, gene %in% liver.genes), aes(y = exprs, x = tissgeno)) + geom_boxplot()

dat.delta <- subset(dat.mean, geno == "SV129") %>%
  group_by(geno, gene) %>%
  summarise(exprs.delta = exprs[1] - exprs[2])  # liver - kidney

ggplot(subset(dat.delta, gene %in% liver.genes), aes(x = exprs.delta)) + geom_histogram(bins = 40) + theme_bw()


# Load MARA results -------------------------------------------------------

jmain <- "/home/yeung/data/tissue_specificity/mara_results"
outmain <- paste0(jmain, "/mara_outputs", suffix)
outdir <- file.path(outmain, paste0("center.", do.center, suffix))
maraoutdir <- file.path(outdir, E.subdir)

act.s <- LoadActivitiesLong(indir = maraoutdir, shorten.motif.name = TRUE, make.cnames = FALSE)
act.s <- MakeCnamesLivKidWTKO(act.s)

fourier.scale <- 4
zscore.min <- 1.25
omega <- 2 * pi / 24
hr.shift <- 3
mrna.hl <- log(2) / (omega / tan(omega * hr.shift))  # convert hour shift to half-life 
comp <- 1

act.s.complex <- ProjectWithZscore(act.s, omega, n = fourier.scale)
sig.motifs <- unique(as.character(subset(act.s.complex, zscore > zscore.min)$gene))
s.act <- SvdOnComplex(subset(act.s.complex, gene %in% sig.motifs), value.var = "exprs.transformed")

# plot hits
act.s.shift <- act.s
act.s.shift$time <- act.s$time - hr.shift
act.s.shift$time <- sapply(act.s.shift$time, function(x) ifelse(x < 0, x + 48, x))

eigens.act <- GetEigens(s.act,
                        eigenval = FALSE,
                        period = 24, 
                        comp = comp, 
                        adj.mag = TRUE, 
                        constant.amp = 6, 
                        label.n = 20, 
                        jtitle = "", 
                        peak.to.trough = TRUE, 
                        dot.col = "black", 
                        dotsize = 2, 
                        dotshape = 18, 
                        label.gene = c("ELF1.2.4"),
                        half.life = mrna.hl)



# Plot top hits -----------------------------------------------------------

# # targets of ROR
# jmotif <- "RORA"


motif.zscore <- hash(as.character(subset(act.s.complex, tissue == jtissgeno)$gene), signif(subset(act.s.complex, tissue == jtissgeno)$zscore, digits = 2))
tfs <- GetTFs(get.mat.only = TRUE, shorten.motif.name = TRUE)
# rownames(tfs) <- sapply(rownames(tfs), RemoveP2)


# Plot tons of things -----------------------------------------------------


# Do AR separately (include target.genes)
plotmain <- "/home/yeung/projects/tissue-specificity/plots/mara_liver_kidney_modules_on_liverDHS"
plotf <- file.path(plotmain, paste0("target_gene_analysis", suffix, ".", jmotif, ".targetgenes.pdf"))
pdf(plotf)

jzscore <- motif.zscore[[jmotif]]
jmotif.title <- paste(jmotif, "Zscore:", jzscore, sep = " ")
print(PlotActivitiesWithSE(subset(act.s.shift, gene == jmotif), jtitle = jmotif.title) + theme_bw() + theme(aspect.ratio = 1))

genes.all <- unlist(sapply(jmotif, GetGenesFromMotifs, tfs))
genes.that.fit <- genes.all  # take all

for (gene.hit in genes.that.fit){
  if (jmod %in% c("Liver_SV129,Liver_BmalKO", "Liver_SV129", "Liver_BmalKO")){
    gene.prot <- gene.hit
    jprot.long <- prot.long
  } else {
    gene.prot <- ""
    jprot.long <- NA
  }
  print(PlotGeneTissuesWTKO(subset(dat.wtko, gene == gene.hit), jtitle = gene.hit))
  print(PlotmRNAActivityProtein(dat.wtko, act.s.shift, gene.dat = gene.hit, prot.long = jprot.long, gene.act = jmotif, gene.prot = gene.prot, jtiss = jtiss.nogeno, dotsize = 3, themesize = 22, single.day = TRUE) + theme(strip.text = element_blank()))
}
Nsub <- subset(N.sum, motif == jmotif) %>% 
  arrange(desc(sitecount))
print(ggplot(Nsub, aes(x = log2(sitecount))) + geom_histogram(bins = 20) + ggtitle(paste("Sitecounts of", jmotif.title)) + theme_bw())
target.genes <- Nsub$gene[1:top.n]

sink(file = file.path(plotmain, paste0(jmotif, suffix, ".target_genes.txt")))
for (g in target.genes){
  cat(g); cat("\t"); cat(phases[[g]]); cat("\n")
}
sink()
print(ggplot(subset(Nsub, gene %in% target.genes), aes(x = mean)) + geom_histogram(bins = 25) + ggtitle(paste("Means Motif:", jmotif.title)) + theme_bw())
circular_phase24H_histogram(subset(Nsub, gene %in% target.genes)$phase, jtitle = paste("Phase Motif", jmotif.title), cex.axis = jcex, cex.lab = jcex, color_hist = "red")

# plot target genes
for (g in target.genes){
  print(PlotGeneTissuesWTKO(subset(dat.wtko, gene == g), jtitle = g))
}

dev.off()


# plotmain <- "/home/yeung/projects/tissue-specificity/plots/mara_liver_kidney_modules_on_liverDHS"
# plotf <- file.path(plotmain, paste0("target_gene_analysis", suffix, ".pdf"))
# pdf(plotf, useDingbats = TRUE)
# 
# print(eigens.act$u)
# print(eigens.act$v)
# 
# for (jmotif in sig.motifs){
#   jzscore <- motif.zscore[[jmotif]]
#   jmotif.title <- paste(jmotif, "Zscore:", jzscore, sep = " ")
#   print(PlotActivitiesWithSE(subset(act.s.shift, gene == jmotif), jtitle = jmotif.title) + theme_bw() + theme(aspect.ratio = 1))
#   
#   genes.all <- unlist(sapply(jmotif, GetGenesFromMotifs, tfs))
#   genes.that.fit <- genes.all  # take all
#   
#   for (gene.hit in genes.that.fit){
#     if (jmod %in% c("Liver_SV129,Liver_BmalKO", "Liver_SV129", "Liver_BmalKO")){
#       gene.prot <- gene.hit
#       jprot.long <- prot.long
#     } else {
#       gene.prot <- ""
#       jprot.long <- NA
#     }
#     print(PlotGeneTissuesWTKO(subset(dat.wtko, gene == gene.hit), jtitle = gene.hit))
#     print(PlotmRNAActivityProtein(dat.wtko, act.s.shift, gene.dat = gene.hit, prot.long = jprot.long, gene.act = jmotif, gene.prot = gene.prot, jtiss = jtiss.nogeno, dotsize = 3, themesize = 22, single.day = TRUE) + theme(strip.text = element_blank()))
#   }
#   Nsub <- subset(N.sum, motif == jmotif) %>% 
#     arrange(desc(sitecount))
#   print(ggplot(Nsub, aes(x = log2(sitecount))) + geom_histogram(bins = 20) + ggtitle(paste("Sitecounts of", jmotif.title)) + theme_bw())
#   target.genes <- Nsub$gene[1:top.n]
#   print(ggplot(subset(Nsub, gene %in% target.genes), aes(x = mean)) + geom_histogram(bins = 25) + ggtitle(paste("Means Motif:", jmotif.title)) + theme_bw())
#   circular_phase24H_histogram(subset(Nsub, gene %in% target.genes)$phase, jtitle = paste("Phase Motif", jmotif.title), cex.axis = jcex, cex.lab = jcex, color_hist = "red")
# }
# dev.off()


# 
# 
# # targets of Ebox
# jmotif <- "bHLH_family.p2"
# Nsub <- subset(N.sum, motif == jmotif) %>% 
#   arrange(desc(sitecount))
# ggplot(Nsub, aes(x = log2(sitecount))) + geom_histogram(bins = 20) + ggtitle(paste("Sitecounts of", jmotif))
# target.genes <- Nsub$gene[1:top.n]
# ggplot(subset(Nsub, gene %in% target.genes), aes(x = mean)) + geom_histogram(bins = 25) + ggtitle(paste("Means Motif:", jmotif))
# circular_phase24H_histogram(subset(Nsub, gene %in% target.genes)$phase, jtitle = paste("Phase Motif", jmotif), cex.axis = jcex, cex.lab = jcex, color_hist = "red")
# 
# # targets of AR
# jmotif <- "AR.p2"
# Nsub <- subset(N.sum, motif == jmotif) %>% 
#   arrange(desc(sitecount))
# ggplot(Nsub, aes(x = log2(sitecount))) + geom_histogram(bins = 20) + ggtitle(paste("Sitecounts of", jmotif))
# target.genes <- Nsub$gene[1:top.n]
# ggplot(subset(Nsub, gene %in% target.genes), aes(x = mean)) + geom_histogram(bins = 25) + ggtitle(paste("Means Motif:", jmotif))
# circular_phase24H_histogram(subset(Nsub, gene %in% target.genes)$phase, jtitle = paste("Phase Motif", jmotif), cex.axis = jcex, cex.lab = jcex, color_hist = "red")
# 
# # targets of FOXA2
# jmotif <- "FOXA2.p3"
# Nsub <- subset(N.sum, motif == jmotif) %>% 
#   arrange(desc(sitecount))
# ggplot(Nsub, aes(x = log2(sitecount))) + geom_histogram(bins = 20) + ggtitle(paste("Sitecounts of", jmotif))
# target.genes <- Nsub$gene[1:top.n]
# ggplot(subset(Nsub, gene %in% target.genes), aes(x = mean)) + geom_histogram(bins = 25) + ggtitle(paste("Means Motif:", jmotif))
# circular_phase24H_histogram(subset(Nsub, gene %in% target.genes)$phase, jtitle = paste("Phase Motif", jmotif), cex.axis = jcex, cex.lab = jcex, color_hist = "red")
