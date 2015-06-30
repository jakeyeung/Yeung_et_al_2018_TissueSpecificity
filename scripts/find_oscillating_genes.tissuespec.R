# 2015-06-18
# find_oscillating_genes.kallistoarray.R
# We want to find oscillating genes, with p-value and amplitude proportional 
# to a core clock gene as a cutoff.

library(ggplot2)
library(dplyr)

setwd("/home/yeung/projects/tissue-specificity")

# Functions ---------------------------------------------------------------

source("scripts/functions/PlotGeneTpm.R")
source("scripts/functions/LoadArrayRnaSeq.R")
source("scripts/functions/LoadKallisto.R")
source("scripts/functions/GetTissueTimes.R")
source("scripts/functions/FitRhythmic.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/WriteGeneListMat.R")
source("scripts/functions/TissueSpecificRhythmicsFunctions.R")

# Load matrix -------------------------------------------------------------

# array.path <- "data/exprs_matrices/array_adj_to_kallisto.slope07.txt"
# rna.seq.path <- "/home/yeung/projects/tissue-specificity/data/kallisto/abundance.genecounts.matrix.txt"
# ka.long <- LoadLong(array.path, rna.seq.path, scale.factor = 100, pseudocount = 1)
ka.long <- LoadArrayRnaSeq()
tpm.long <- LoadKallisto()

PlotGeneAcrossTissues(subset(ka.long, gene == "Arntl"))


# Fit rhythmic ------------------------------------------------------------

fits <- FitRhythmicDatLong(ka.long)

fits.relamp <- GetRelamp(fits = fits, max.pval = 1e-3)

# Find rhythmic genes ----------------------------------------------

pval.min <- 1e-5
pval.max <- 0.05
relamp.max <- 0.1
mean.cutoff <- 4

fits.relamp$is.rhythmic <- apply(as.matrix(fits.relamp), 1, IsRhythmicApply, 
                          pval.min = pval.min, pval.max = pval.max, relamp.max = relamp.max, cutoff = mean.cutoff)


# Do not consider Cere, Hypo and BS into this -----------------------------

filter.tissues <- c("Cere", "BS", "Hypo", "Adr", "Aorta", "BFAT", "Heart", "Lung", "Mus", "WFAT")

fits.relamp.filt <- subset(fits.relamp, !tissue %in% filter.tissues)


# Get tissue-specific genes -----------------------------------------------

# only tissue-specific if it contains TRUE and FALSE in a gene
fits.relamp.filt <- fits.relamp.filt %>%
  group_by(gene) %>%
  do(IsTissueSpecificDplyr(.))

tissue.spec <- unique(subset(fits.relamp.filt, is.tiss.spec == TRUE)$gene)


# Sort by significance ----------------------------------------------------

fits.tiss.spec <- subset(fits.relamp.filt, is.tiss.spec == TRUE)
head(data.frame(fits.tiss.spec[order(fits.tiss.spec$pval), ]), n = 50)
head(data.frame(fits.tiss.spec[order(fits.tiss.spec$relamp, decreasing = TRUE), ]), n = 50)

tissue.spec.amp <- unique(data.frame(fits.tiss.spec[order(fits.tiss.spec$relamp, decreasing = TRUE), ])$gene)


# Plot by amplitude -------------------------------------------------------

tpm.filt <- subset(tpm.long, gene_name %in% tissue.spec)
fits.filt <- subset(fits.relamp.filt, gene %in% tissue.spec)
dat.filt <- subset(ka.long, gene %in% tissue.spec)
tpm.genes <- unique(tpm.filt$gene)

outpath <- "plots/tissue_specific_rhythmic_genes/liver_v_kidney.pdf"
print(paste("Plotting liver vs kidney to", outpath))
start.time <- Sys.time()
pdf(outpath)
for (jgene in tissue.spec.amp){
  if (!jgene %in% tpm.genes) next
  print(PlotGeneAcrossTissues(subset(dat.filt, gene == jgene)))
  print(PlotTpmAcrossTissues(subset(tpm.filt, gene_name == jgene), log2.transform = FALSE))
}
dev.off()
print(Sys.time() - start.time)


# Plot examples -----------------------------------------------------------

# by amp
jgene <- "Dnmt3b"
jgene <- "BC029214"
jgene <- "Lrfn3"
jgene <- "Atf5"
jgene <- "Gfod2"
jgene <- "Mafb"
jgene <- "Ddit4"
jgene <- "Gys2"  # good dist
jgene <- "Eps8l2"  # hnf regulated? Adr spikey. Almost identical DHS sites btw kid and liv
jgene <- "Slc25a47"
jgene <- "Fads3"
jgene <- "Mapk15"
jgene <- "Pla2g12a"
jgene <- "Tfec"
jgene <- "Pcp4l1"  # OK distance
jgene <- "Ggt6"
jgene <- "Maff"
jgene <- "Homer2"  # distance good but lowly expressed
jgene <- "Crip2"
jgene <- "Zc3h12d"  # possible candidate
jgene <- "Flcn"
jgene <- "Rbm12b1"
jgene <- "Pxmp4"
jgene <- "Uck2"
jgene <- "Rcan1"  # good candidate for alternative promoter
jgene <- "Rnd1"
jgene <- "Gpr97"
jgene <- "Slc38a2"  # good distance with eRNA evidence
jgene <- "Neu3"  # good distance with ERNA evidence
jgene <- "Spcs2"  # neu3's neighbour
jgene <- "Rad9a"  
jgene <- "Ppp1ca"  # Rad9a's neighbhour
jgene <- "Pik3ap1"  # good candidate
jgene <- "Hp"
jgene <- "Steap3"
jgene <- "Cep44"

# by pval
jgene <- "Slc44a1"  # good distance + eRNA evidence
jgene <- "Hsd3b7"
jgene <- "Rbm12b1"
jgene <- "Sco2"
jgene <- "Pdgfc"  # good distance, no eRNA evidence
jgene <- "Tjp3"
jgene <- "Repin1"
jgene <- "Fbxo6"
jgene <- "Slc4a4"  # good distance, alt promoter, has eRNA
jgene <- "Lgals9"
jgene <- "Nhlrc2"  # the double promoter
jgene <- "Lrp5"  # Adr + Liver rhythmic, maybe a candidate but it's muddy


PlotGeneAcrossTissues(subset(ka.long, gene == jgene))
PlotTpmAcrossTissues(subset(tpm.filt, gene_name == jgene), log2.transform = TRUE)


# Plot candidates ---------------------------------------------------------

candidates <- 
  c("Slc44a1",
"Pik3ap1",
"Neu3",
"Slc38a2",
"Rcan1",
"Slc4a4",
"Gys2",
"Homer2",
"Zc3h12d",
"Pdgfc")

pdf("plots/tissue_specific_rhythmic_genes/candidates_4c_for_jerome.pdf")
for (g in candidates){
  print(PlotGeneAcrossTissues(subset(ka.long, gene == g)))
  print(PlotTpmAcrossTissues(subset(tpm.filt, gene_name == g & tissue %in% c("Liver", "Kidney")), log2.transform = FALSE))
}
dev.off()
