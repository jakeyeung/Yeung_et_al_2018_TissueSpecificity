# 2016-09-02
# Jake Yeung
# See how activity changes over different lambdas


rm(list=ls())

library(dplyr)
library(ggplot2)
library(hash)

source("scripts/functions/ListFunctions.R")
source("scripts/functions/LiverKidneyFunctions.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/NcondsFunctions.R")
source("scripts/functions/SvdFunctions.R")
source("scripts/functions/GetClockGenes.R")
source("scripts/functions/BiomartFunctions.R")
source("scripts/functions/PlotActivitiesFunctions.R")
source("scripts/functions/GetTFs.R")


# Functions ---------------------------------------------------------------

PlotComplex.coloredbylambda <- function(dat, omega = 2 * pi / 24, 
                           title = "My title", xlab = "Amplitude of activity", ylab = "Phase of activity (CT)", 
                           ampscale = 2, constant.amp = 5, dot.col = "gray85", jsize = 22){
    # Convert complex to amplitude (2 * fourier amplitude) and phase, given omega.
    # then plot in polar coordinates
    # fourier amplitudes are half-amplitudes of the sine-wave
    # http://www.prosig.com/signal-processing/FourierAnalysis.pdf
    # Default: ampscale = 2 because fourier amplitude is more like half-amplitude by default
    # 
    # Add colours
    # http://stackoverflow.com/questions/20808009/remove-extra-space-and-ring-at-the-edge-of-a-polar-plot
  
    # dat <- data.frame(amp = Mod(vec.complex) * ampscale,
    #                  phase = ConvertArgToPhase(Arg(vec.complex), omega = omega),
    #                  label = labels)
    
    amp.max <- ceiling(max(dat$amp) * 2) / 2
    if (amp.max <= 1){
      amp.step <- 0.5
    } else {
      amp.step <- 1
    }
    m <- ggplot(data = dat, aes(x = amp, y = phase, label = gene, colour = lambda)) + 
      geom_point(size = 0.5, colour = dot.col) +
      coord_polar(theta = "y") + 
      xlab(xlab) +
      ylab(ylab) +
      ggtitle(title) +
      scale_y_continuous(limits = c(0, 24), breaks = seq(6, 24, 6)) + 
      scale_x_continuous(limits = c(0, amp.max), breaks = seq(0, amp.max, length.out = 2)) + 
      theme_bw(jsize) + 
      geom_vline(xintercept = seq(0, amp.max, length.out = 2), colour = "grey50", size = 0.2, linetype = "dashed") +
      geom_hline(yintercept = seq(6, 24, by = 6), colour = "grey50", size = 0.2, linetype = "solid") +
      theme(panel.grid.major = element_line(size = 0.5, colour = "grey"), panel.grid.minor = element_blank(), 
            panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position="bottom",
            panel.border = element_blank(),
            legend.key = element_blank(),
            axis.ticks = element_blank(),
            panel.grid  = element_blank())
    
    # add text
    dat.txt <- subset(dat, gene != "")
    if (constant.amp != FALSE){
      m <- m + geom_text_repel(data = dat.txt, aes(x = amp, y = phase, label = gene), size = constant.amp)
    } else {
      m <- m + geom_text_repel(data = dat.txt, aes(x = amp, y = phase, size = amp, label = gene))
    }
    return(m)
}

# Load --------------------------------------------------------------------

outmain.lambda <- "/home/yeung/projects/tissue-specificity/results/MARA.liver_kidney.globalLambda"
outdirs <- list.dirs(outmain.lambda, recursive = FALSE, full.names = TRUE)

outdir <- file.path(outmain.lambda, "promoters.Liver_SV129,Kidney_SV129.Liver_BmalKO,Kidney_BmalKO-Liver_SV129,Kidney_SV129.g=1001.lambda.0.0479929")
omega <- 2 * pi / 24

act.complex <- lapply(outdirs, function(outdir){
  indir <-  file.path(outdir, "atger_with_kidney.bugfixed")
  source("scripts/functions/LoadActivitiesLong.R")
  act.long <- LoadActivitiesLongKidneyLiver(indir, collapse.geno.tissue = TRUE, shorten.motif.name = TRUE)
  # use Liver_SV129, but should be global lambda so it shouldnt matter
  lambda.f <- file.path(indir, "Liver_SV129", "Lambda")
  lambda <- as.numeric(unlist(read.table(lambda.f)))
  act.complex <- act.long %>%
    group_by(gene, tissue) %>%
    do(ProjectToFrequency2(., omega, add.tissue=TRUE))
  act.complex <- act.complex %>%
    mutate(amp = Mod(exprs.transformed) * 4) %>%
    mutate(phase = ConvertArgToPhase(Arg(exprs.transformed), omega = omega))
  act.complex$lambda <- lambda  
  return(act.complex)
})
act.complex <- do.call(rbind, act.complex)


# Show how HIC1 amplitude changes over time -------------------------------

jsub <- subset(act.complex, tissue == "Liver_SV129" & gene %in% c("HIC1", "bHLH_family", "RORA", "NFIL3"))


ggplot(jsub, aes(x = lambda, y = amp)) + geom_point()  + facet_wrap(~gene)

PlotComplex.coloredbylambda(jsub, xlab = "Log2 Fold Change", ylab = "ZT", title = "Global lambda", jsize = 16, constant.amp = 4)

