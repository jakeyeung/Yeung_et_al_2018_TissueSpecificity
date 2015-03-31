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

ProjectToFrequency2 <- function(df, omega){
  # simpler than ProjectToFrequency().
  # expect df to be gene i in tissue c with column names time and exprs
  exprs.transformed <- DoFourier(df$exprs, df$time, omega = omega)
  return(data.frame(exprs.transformed = exprs.transformed))
}

ConvertArgToPhase <- function(phase.rads, omega){
  # convert phase in radians to phase in time, using omega.
  # expect phase to be between -pi to pi, convert that
  # to 0 to 24 hours.
  
  # convert range from -pi to pi to 0 to 2pi
  phase.rads[which(phase.rads < 0)] <- phase.rads[which(phase.rads < 0)] + 2 * pi
  phase <- phase.rads / omega
  return(phase)
}

PlotComplex2 <- function(vec.complex, labels, omega = 2 * pi / 24, title = "My title"){
  # Convert complex to amplitude (2 * fourier amplitude) and phase, given omega.
  # then plot in polar coordinates
  # fourier amplitudes are half-amplitudes of the sine-wave
  # http://www.prosig.com/signal-processing/FourierAnalysis.pdf
  df <- data.frame(amp = Mod(vec.complex) * 2,
                   phase = ConvertArgToPhase(Arg(vec.complex), omega = omega),
                   label = labels)
  m <- ggplot(data = df, aes(x = amp, y = phase, label = label)) + 
        geom_point(size = 0.5) +
        coord_polar(theta = "y") + 
        xlab("Amplitude of activity") +
        ylab("Phase of activity") +
        geom_text(aes(x = amp, y = phase, size = amp), hjust = 0, vjust = 0) +
        ggtitle(title) +
      scale_y_continuous(limits = c(0, 24), breaks = seq(2, 24, 2))
}

PlotEigengene <- function(svd.obj, comp, omega = 2 * pi / 24){
  eigengene <- svd.obj$v[, comp]
  v.plot <- PlotComplex2(eigengene, labels = rownames(svd.obj$v), omega = omega, title = paste("Right singular value", comp))
  print(v.plot)
}

PlotEigensamp <- function(svd.obj, comp, omega = 2 * pi / 24){
  eigensamp <- svd.obj$u[, comp]
  u.plot <- PlotComplex2(eigensamp, labels = rownames(svd.obj$u), omega = omega, title = paste("Left singular value", comp))   
  print(u.plot)
}


# Define constants --------------------------------------------------------

remove.wfat <- TRUE

# Load file ---------------------------------------------------------------


# out.fpath from cbind_activities.R
in.fpath <- "/home/yeung/projects/tissue-specificity/results/outputs_all_genes_MARA/MARA.33motifs.corrected/activities.all"

act.all <- read.table(in.fpath)
motif_activities.plotout <- "plots/all_genes.MARA.33motifs.corrected.all_motifs.activities.pdf"


# Create long dataframe ---------------------------------------------------

tissues <- GetTissues(colnames(act.all))
times <- GetTimes(colnames(act.all))
act.long <- data.frame(gene = rep(rownames(act.all), ncol(act.all)),  # i know it's a motif, bare with me.
                       tissue = rep(tissues, each = nrow(act.all) * length(times)),
                       time = as.numeric(rep(times, length(tissues), each = nrow(act.all))),
                       exprs = as.numeric(unlist(act.all)))


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
  
  v.plot <- PlotComplex2(eigengene, labels = colnames(act.proj.wide), omega = omega, title = paste("Right singular value", comp))
  u.plot <- PlotComplex2(eigensamp, labels = rownames(act.proj.wide), omega = omega, title = paste("Left singular value", comp)) 
  print(v.plot)
  print(u.plot)
}
dev.off()
