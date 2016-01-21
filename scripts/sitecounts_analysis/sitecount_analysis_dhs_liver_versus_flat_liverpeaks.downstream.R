# 2015-12-17
# Downstream analysis of comparing DHS sites liver versus flat on liverpeaks

source("scripts/functions/FisherTestSitecounts.R")


# Load --------------------------------------------------------------------

load("Robjs/S.collapse.flat.dist1000.Robj", verbose=T)
load("Robjs/S.collapse.liver.dist1000.Robj", verbose=T)
load("Robjs/N.long.flatliver_genes.all_motifs.dist1000.Robj", verbose=T)


# Do ----------------------------------------------------------------------


# label different models
S.collapse.flat$model <- "Flat"
S.collapse$model <- "Liver"

# merge together
S.collapse <- rbind(S.collapse, S.collapse.flat)

# filter for liver peaks
S.collapse <- subset(S.collapse, peak.type == "Liver")

# load up sitecounts for liver and flat genes

N.long.all <- subset(N.long.all, peak %in% S.collapse$peak)

tests.motif.livflat <- N.long.all  %>%
  group_by(motif) %>%
  do(RunFisherDHS(., S.collapse, jmodel.col = "model"))

# plot volcano

RunFisherDHS(subset(N.long.all, motif == "HIF1A.p2"), S.collapse, jmodel.col = "model", jshow.table = TRUE)

ggplot(tests.motif.livflat, aes(x = odds.ratio, y = -log10(p.value), label = motif)) + geom_text()



# Redo with Robjs ---------------------------------------------------------

load("Robjs/S.collapse.liverflat.liver_peaks_only.Robj", verbose=T)
load("Robjs/N.long.flatliver_genes.all_motifs.dist1000.liver_peaks_only.Robj", verbose=T)

tests.motif.livflat <- N.long.all.livpeaks %>%
  group_by(motif) %>%
  do(RunFisherDHS(., S.collapse.livflat, jmodel.col = "model"))
ggplot(tests.motif.livflat, aes(x = odds.ratio, y = -log10(p.value), label = motif)) + geom_text()


