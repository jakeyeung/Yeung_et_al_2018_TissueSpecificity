# 2017-04-28
# Jake Yeung
# Can we identify rhythmic regulators by filtering liver DHSs then running MARA?

rm(list=ls())

library(ggplot2)
library(reshape2)
library(dplyr)


source("scripts/functions/LiverKidneyFunctions.R")
source("scripts/functions/LoadActivitiesLong.R")
source("scripts/functions/PlotActivitiesFunctions.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/CosSineFunctions.R")
source("scripts/functions/HalfLifeFunctions.R")
source("scripts/functions/SvdFunctions.R")
source("scripts/functions/GetTFs.R")


# Load liver DHS ----------------------------------------------------------

jweight <- 0.8  # take all liver DHSs to all genes in Liver_SV129
jweight <- 0  # take all liver DHSs to all genes in Liver_SV129

flatampmax <- 0.1

promoters.only <- FALSE
all.genes <- FALSE
suffix <- paste0(".weight.", jweight, ".flatampmax.", flatampmax, ".promoters.", promoters.only, ".all_genes.", all.genes)

inf <- paste0("/home/yeung/projects/tissue-specificity/Robjs/three_way_cooccurence/three.way.cooccurrence.bugfixed.nmodels.2.K.200.weight.", jweight, ".MergePeaks.FALSE.nullmodel.JI.flatampmax.", flatampmax, ".withNmatallNmatfreqs.RemoveZeroCounts.Robj")
load(inf, v=T)

jmeth <- "g=1001"
load("Robjs/liver_kidney_atger_nestle/fits.long.multimethod.filtbest.staggeredtimepts.bugfixed.annotated.Robj", v=T)
fits.long.filt <- subset(fits.long.filt, method == jmeth)

liver.genes.all <- as.character(subset(fits.long.filt, model == "Liver_SV129")$gene)

if (all.genes){
  liver.genes <- liver.genes.all
} else {
  liver.peaks <- unique(as.character(subset(N.mat.all, model %in% "rhyth")$peak))
  liver.genes <- unique(as.character(subset(N.mat.all, model %in% "rhyth")$gene))
  print(paste("N peaks:, ", length(liver.peaks)))
  print(paste("N genes:, ", length(liver.genes)))
}


# Get gene expression over time and genotypes -----------------------------

load("Robjs/liver_kidney_atger_nestle/dat.long.liverkidneyWTKO.bugfixed.Robj", v=T); dat.wtko <- dat.long; rm(dat.long)
dat.wtko <- StaggeredTimepointsLivKid(dat.wtko)

# Load MARA output --------------------------------------------------------

maraoutdir <- paste0("/home/yeung/data/tissue_specificity/mara_results/mara_outputs", suffix, "/center.TRUE", suffix, "/centered.TRUE")
act.s <- LoadActivitiesLong(indir = maraoutdir, shorten.motif.name = TRUE, make.cnames = FALSE)
act.s$sampname <- act.s$tissue
act.s$tissue <- as.character(sapply(as.character(act.s$sampname), function(s) strsplit(s, "_")[[1]][[1]]))
act.s$time <- as.numeric(sapply(as.character(act.s$sampname), function(s) strsplit(s, "_")[[1]][[2]]))
act.s$geno <- as.character(sapply(as.character(act.s$sampname), function(s) strsplit(s, "_")[[1]][c(-1, -2)]))
act.s$tissue <- paste(act.s$tissue, act.s$geno, sep = "_")
act.s$tissue <- factor(act.s$tissue, levels = c("Liver_SV129", "Liver_BmalKO", "Kidney_SV129", "Kidney_BmalKO"))
act.s$experiment <- "rnaseq"
act.s$sampname <- NULL

fourier.scale <- 4
zscore.min <- 1.25
omega <- 2 * pi / 24
hr.shift <- 3
mrna.hl <- log(2) / (omega / tan(omega * hr.shift))  # convert hour shift to half-life 
comp <- 1

act.s.complex <- ProjectWithZscore(act.s, omega, n = fourier.scale)
sig.motifs <- unique(as.character(subset(act.s.complex, zscore > zscore.min)$gene))
s.act <- SvdOnComplex(subset(act.s.complex, gene %in% sig.motifs), value.var = "exprs.transformed")

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

eigens.act.fancy.LivWTKO <- GetEigens(s.act, period = 24, comp = comp, adj.mag = TRUE, 
                                      eigenval = FALSE,
                                      constant.amp = 5, 
                                      label.n = Inf, jtitle = "", 
                                      peak.to.trough = TRUE, 
                                      dot.col = "black", 
                                      # dot.col = col.hash, 
                                      dotsize = 6, 
                                      dotshape = 18,
                                      disable.text = FALSE, 
                                      add.arrow = TRUE,
                                      disable.repel = TRUE,
                                      half.life = mrna.hl)

print(eigens.act$v.plot)
print(eigens.act$u.plot)

print(eigens.act.fancy.LivWTKO$v.plot)
print(eigens.act.fancy.LivWTKO$u.plot)

# also plot TF activity adjusted by half-life
s.ampphase <- GetAmpPhaseFromActivities(act.s, mrna.hl, jtiss = jtiss, jgeno = "SV129")
# just shift by 3 hours
act.s.shift <- act.s
act.s.shift$time <- act.s$time - hr.shift
act.s.shift$time <- sapply(act.s.shift$time, function(x) ifelse(x < 0, x + 48, x))

jmotif <- "RORA"
print(PlotActivitiesWithSE(subset(act.s.shift, gene == jmotif), jtitle = jmotif) + theme_bw())


source("scripts/functions/PlotGeneAcrossTissues.R")
load("Robjs/liver_kidney_atger_nestle/fits.bytiss.bugfixed.Robj", v=T)

liv.rhyth <- as.character(subset(fits.bytiss, tissue == "Liver_SV129" & amp > 0.25 & pval < 1e-3)$gene)
tfs <- GetTFs(get.mat.only = TRUE)

# Print genes that match model
jmotifs <- names(head(eigens.act$eigensamp[order(abs(eigens.act$eigensamp), decreasing = TRUE)], n = 20))

for (jmotif in jmotifs){
  genes.all <- unlist(sapply(jmotif, GetGenesFromMotifs, tfs))
  genes.that.fit <- genes.all[which(genes.all %in% liv.rhyth)]
  if (jmotif == "RORA"){
    genes.that.fit <- c(genes.that.fit, "Nr1d1")
  }
  print(paste("Rhyth genes:", genes.that.fit))
  if (length(genes.that.fit) > 0){
    for (gene.hit in genes.that.fit){
      # if (jmod %in% c("Liver_SV129,Liver_BmalKO", "Liver_SV129", "Liver_BmalKO")){
      #   gene.prot <- gene.hit
      #   jprot.long <- prot.long
      # } else {
      #   gene.prot <- ""
      #   jprot.long <- NA
      # }
      print(PlotGeneTissuesWTKO(subset(dat.wtko, gene == gene.hit), jtitle = gene.hit))
      print(PlotActivitiesWithSE(subset(act.s.shift, gene == jmotif), jtitle = jmotif) + theme_bw())
      # print(PlotmRNAActivityProtein(dat.wtko, act.s, gene.dat = gene.hit, prot.long = jprot.long, gene.act = jmotif, gene.prot = gene.prot, jtiss = jtiss, dotsize = 3, themesize = 22) + theme(strip.text = element_blank()))
      # print(PlotmRNAActivityProtein(dat.wtko, act.s.shift, gene.dat = gene.hit, prot.long = jprot.long, gene.act = jmotif, gene.prot = gene.prot, jtiss = jtiss, dotsize = 3, themesize = 22) + theme(strip.text = element_blank()))
      # print(PlotmRNAActivityProtein(dat.wtko, act.s.shift, gene.dat = gene.hit, prot.long = jprot.long, gene.act = jmotif, gene.prot = gene.prot, jtiss = jtiss, dotsize = 3, themesize = 22, single.day = TRUE) + theme(strip.text = element_blank()))
      # print(PlotmRNAActivityProtein(dat.wtko, s.ampphase, gene.dat = gene.hit, prot.long = jprot.long, gene.act = jmotif, gene.prot = gene.prot, jtiss = jtiss, dotsize = 3, themesize = 22, act.in.sine.cos = TRUE) + theme(strip.text = element_blank()))
      # print(PlotmRNAActivityProtein(dat.wtko, s.ampphase, gene.dat = gene.hit, prot.long = jprot.long, gene.act = jmotif, gene.prot = gene.prot, jtiss = jtiss, dotsize = 3, themesize = 22, act.in.sine.cos = TRUE, single.day = TRUE) + theme(strip.text = element_blank()))
      # print(PlotmRNAActivityProtein(dat.wtko, act.s, gene.dat = gene.hit, prot.long = jprot.long, gene.act = jmotif, gene.prot = gene.prot, jtiss = "both", dotsize = 2, themesize = 14) + theme(strip.text = element_blank()))
    }
  }
}

# 
# # Old ---------------------------------------------------------------------
# 
# 
# mara <- LoadMaraOutput(maraoutdir)
# 
# zscores <- mara$zscores
# act.long <- mara$act.long
# act.long$samp <- NULL
# 
# act.long$tissue <- sapply(as.character(act.long$sample), function(s) strsplit(s, "_")[[1]][[1]])
# act.long$time <- as.numeric(sapply(as.character(act.long$sample), function(s) strsplit(s, "_")[[1]][[2]]))
# act.long$geno <- as.character(sapply(as.character(act.long$sample), function(s) strsplit(s, "_")[[1]][[3]]))
# 
# PlotGeneTissuesWTKO(subset(act.long, gene == "RORA.p2"))
# PlotGeneTissuesWTKO(subset(act.long, gene == "bHLH_family.p2"))
# PlotGeneTissuesWTKO(subset(act.long, gene == "TFAP4.p2"))
# PlotGeneTissuesWTKO(subset(act.long, gene == "HNF1A.p2"))
# PlotGeneTissuesWTKO(subset(act.long, gene == "FOSL2.p2"))
# PlotGeneTissuesWTKO(subset(act.long, gene == "AR.p2"))
# PlotGeneTissuesWTKO(subset(act.long, gene == "FOXA2.p3"))
# 
# 
