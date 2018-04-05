# plot_activities_across_tissues.R
# Jake Yeung
# 2015-03-30
# After running cbind_activities.R, we get a matrix of activities across all 
# tissues. Plot polar coordinate figure of the activities. Possibly even do SVD analysis

library(plyr)

# Functions ---------------------------------------------------------------

source("scripts/functions/GetTissueTimes.R")
source("scripts/functions/SvdFunctions.R")
source("scripts/functions/ConvertLongToWide.R")
source("scripts/functions/PlotActivitiesFunctions.R")

# Define constants --------------------------------------------------------

remove.wfat <- TRUE

# Load file ---------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
# out.fpath from cbind_activities.R
in.fpath <- "/home/yeung/projects/tissue-specificity/results/MARA/MARA_N_centered_with_SE/multiple/expressed_genes_threshold5/activities/activities.all"
err.fpath <- "/home/yeung/projects/tissue-specificity/results/MARA/MARA_N_centered_with_SE/multiple/expressed_genes_threshold5/activities/standarderrors.all"
# in.fpath <- args[1]


# motif_activities.plotout <- "plots/all_genes.MARA_with_SE.pdf"
# motif_activities.plotout <- args[2]

# Create long dataframe ---------------------------------------------------

act.all <- read.table(in.fpath)
err.all <- read.table(err.fpath)

tissues <- GetTissues(colnames(act.all))
times <- GetTimes(colnames(act.all))
act.long <- data.frame(gene = rep(rownames(act.all), ncol(act.all)),  # i know it's a motif, bare with me.
                       tissue = rep(tissues, each = nrow(act.all) * length(times)),
                       time = as.numeric(rep(times, length(tissues), each = nrow(act.all))),
                       exprs = as.numeric(unlist(act.all)),
                       se = as.numeric(unlist(err.all)))


# Test plot ---------------------------------------------------------------

jgene <- "HNF4A_NR2F1.2.p2"
jgene <- "MEF2.A.B.C.D..p2"
jgene <- "NR1H4.p2"
jgene <- "PITX1..3.p2"
jgene <- "HSF1.2.p2"
jgene <- "bHLH_family.p2"
jgene <- "NFIL3.p2"
jgene <- "DBP.p2"
jgene <- "FOXA2.p3"
jgene <- "GATA6.p2"
jgene <- "REST.p3"
jgene <- "HNF1A.p2"
jgene <- "RORA.p2"
jgene <- "E2F1..5.p2"
jgene <- "NRF1.p2"
jgene <- "ZEB1.p2"
jgene <- "PPARG.p2"
ggplot(data = subset(act.long, gene == jgene) , aes(x = time, y = exprs)) + 
  geom_line() +
  geom_errorbar(aes(ymax = exprs + se, ymin = exprs - se)) +
  geom_errorbar(limits, width=0.25) + 
  facet_wrap(~tissue) + 
  xlab("CT") + 
  ylab("Activity") + 
  ggtitle(jgene)


# Convert to Fourier ------------------------------------------------------

# split
act.split <- split(act.long, act.long$tissue)

# REMOVE WFAT
if (remove.wfat){
  act.split$WFAT <- NULL
}

omega <- 2 * pi / 24
act.split.proj <- lapply(act.split, function(df){
  ddply(df, .(gene), ProjectToFrequency2, omega = omega)
})

# add tissue information
for (tissue in names(act.split.proj)){
  act.split.proj[[tissue]]$tissue <- tissue
}

act.proj <- do.call(rbind, act.split.proj)

act.proj.wide <- ConvertLongToWide(long.df = act.proj, measurement.var = "exprs.transformed")

# Run SVD -----------------------------------------------------------------

act.svd <- svd(act.proj.wide) 

# add row and colnames
rownames(act.svd$u) <- rownames(act.proj.wide)
rownames(act.svd$v) <- colnames(act.proj.wide)



# Plot across tissues -----------------------------------------------------

# plot motif activites across tissues

act.proj$phase <- ConvertArgToPhase(Arg(act.proj$exprs.transformed), omega)
act.proj$amp <- 2 * Mod(act.proj$exprs.transformed)

theme_set(theme_gray(base_size = 14))
pdf(file = motif_activities.plotout)
m <- ggplot(data = act.proj, aes(x = amp, y = phase, label = gene)) + 
      geom_point(size = 0.5) +
      coord_polar(theta = "y") + 
      xlab("Amplitude of activity") +
      ylab("Phase of activity") +
      geom_text(aes(x = amp, y = phase, size = amp), hjust = 0, vjust = 0) +
      scale_y_continuous(limits = c(0, 24), breaks = seq(2, 24, 2)) + 
      scale_size(range=c(0.5, 3)) +
      facet_wrap(~tissue)
print(m)

# Plot first component ----------------------------------------------------


# screeplot
plot(act.svd$d ^ 2, type = 'o')  # eigenvalues

theme_set(theme_gray())
for (comp in seq(1, ncol(act.proj.wide))){
  eigengene <- act.svd$v[, comp]
  eigensamp <- act.svd$u[, comp]
  # rotate to phase of largest magnitude in sample of eigengene
  phase.reference <- Arg(eigengene[which(Mod(eigengene) == max(Mod(eigengene)))])
  rotate.factor <- complex(modulus = 1, argument = phase.reference)
  # rotate eigengene by -phase ref
  eigengene <- eigengene * Conj(rotate.factor)
  # rotate eigensamp by +phase ref
  eigensamp <- eigensamp * rotate.factor
  
  v.plot <- PlotComplex2(eigengene, labels = colnames(act.proj.wide), omega = omega, title = paste("Right singular value", comp))
  u.plot <- PlotComplex2(eigensamp, labels = rownames(act.proj.wide), omega = omega, title = paste("Left singular value", comp)) 
  print(v.plot)
  print(u.plot)
}
dev.off()
