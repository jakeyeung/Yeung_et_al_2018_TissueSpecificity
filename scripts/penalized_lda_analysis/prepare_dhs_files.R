# 2016-02-17
# Prepare DHS and TATA box files


source("scripts/functions/PlotGeneAcrossTissues.R")

# Load --------------------------------------------------------------------

load("Robjs/N.long.flatliver_genes.all_motifs.dist1000.liver_peaks_only.Robj", verbose=T)
load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.15.phase_sd_maxdiff_avg.Robj", verbose=T)
load("Robjs/S.long.Robj", verbose=T)
load("Robjs/S.tissuecutoff.Robj", verbose=T)
load("Robjs/dat.long.fixed_rik_genes.Robj", verbose=T)

# Choose the most rhythmic liver genes to be considered -------------------

liv.genes.all <- unique(N.long.all.livpeaks$gene)

amps <- subset(fits.best, gene %in% liv.genes.all & model == "Liver")$amp.avg 
phases <- subset(fits.best, gene %in% liv.genes.all & model == "Liver")$phase.avg

cutoff <- 0.5  # from density plot
plot(density(amps[!is.na(amps)]))
abline(v = cutoff)

plot(density(phases))

# plot amplitude of to show not much difference btwn early and late phase
plot(density(amps[which(phases < 10)]), xlim = c(0, 1))
plot(density(amps[which(phases > 10)]), xlim = c(0, 1))

liv.genes.top <- subset(fits.best, gene %in% liv.genes.all & model == "Liver" & amp.avg > cutoff)


# Choose flat genes -------------------------------------------------------

flat.genes <- unique(subset(N.long.all.livpeaks, model == "Flat")$gene)

