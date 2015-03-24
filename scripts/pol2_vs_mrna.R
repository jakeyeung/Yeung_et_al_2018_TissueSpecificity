# Compare pol2 activity with mrna.
# pol2_vs_mrna.R
# Jake Yeung
# 2015-03-24

library(devtools)
library(f24.R2.cycling)

# Functions ---------------------------------------------------------------

source("scripts/functions/LoadPol2Activities.R")
source("scripts/functions/LoadArrayRnaSeq.R")
source("scripts/functions/FitRhythmic.R")
source("scripts/functions/PlotFunctions.R")

AvgPhase <- function(phases, amps){
  # weighted average of phase based on its amplitude.
  # larger amp -> larger weight 
  amps.norm <- amps / sum(amps)
  phase.avg <- sum(phases * amps.norm)
  return(phase.avg)
}

AverageAmpPhase <- function(df){
  # Expect amp and phase column names
  df.out <- data.frame(amp = mean(df$amp),
                       phase = AvgPhase(df$phase, df$amp),
                       data = "PolII")
  return(df.out)
}

PlotMrnaVsPol2 <- function(df, textsize=4){
  ggplot(data = df, aes(x = amp, y = phase, label = gene)) + 
    geom_point() + 
    geom_text(size = textsize) + 
    coord_polar(theta = "y") +
    scale_y_continuous(limits = c(0, 24), breaks = seq(2, 24, 2)) + 
    facet_wrap(~data)  # mrna or pol2
}

# Main --------------------------------------------------------------------

cutoff <- 1e-5

pol2.filtered <- LoadPol2Activities()  # $signal and $sitecounts
pol2.signal.all <- LoadPol2Activities.all()
load("/home/yeung/projects/tissue-specificity/results/fits.Robj")  # fit.list from find.oscillating.genes.R to filter rhythmic gene

genes.pol2 <- rownames(pol2.filtered$signal)
pol2.filtered$signal <- data.frame(gene = genes.pol2, pol2.filtered$signal)

dat.long <- LoadArrayRnaSeq()
dat.liver <- subset(dat.long, tissue == "Liver")

genes.all <- unique(dat.liver$gene)

genes.rhythmicl <- sapply(genes.all, function(g){
  pval <- fit.list[["Liver"]][[g]]$pval
  if (is.nan(pval)){
    return(FALSE)
  }
  if (pval < cutoff){
    return(TRUE)
  } else {
    return(FALSE)
  }
})

genes.rhythmic <- genes.all[which(genes.rhythmicl)]

dat.liver.filt <- subset(dat.liver, gene %in% genes.rhythmic)

dat.liver.fit.all <- ddply(dat.liver, .(gene), FitRhythmic, .progress = "text")  # about 3 minutes

dat.liver.filt.fit <- ddply(dat.liver.filt, .(gene), FitRhythmic)


# create long vector
dat.liver.mrna <- data.frame(gene = dat.liver.filt.fit$gene,
                             amp = dat.liver.filt.fit$amp, 
                             phase = dat.liver.filt.fit$phase,
                             data = rep("mRNA", nrow(dat.liver.filt.fit)))
dat.liver.pol2 <- data.frame(gene = pol2.filtered$signal$gene,
                             amp = pol2.filtered$signal$amp,
                             phase = pol2.filtered$signal$phase,
                             data = rep("PolII", nrow(pol2.filtered$signal)))
# take care of duplicates in pol2 data
dat.liver.pol2 <- ddply(dat.liver.pol2, .(gene), AverageAmpPhase)

dat.liver.rhythmic <- rbind(dat.liver.mrna, dat.liver.pol2)
head(dat.liver.rhythmic)

# plot genes rhythmic in pol2
dat.liver.pol2rhythmic <- subset(dat.liver.rhythmic, gene %in% genes.pol2)
# plot genes rhythmic in mrna
dat.liver.mrnarhythmic <- subset(dat.liver.rhythmic, gene %in% genes.rhythmic)

PlotMrnaVsPol2(dat.liver.pol2rhythmic)
PlotMrnaVsPol2(dat.liver.mrnarhythmic)
