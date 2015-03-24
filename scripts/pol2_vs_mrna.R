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
# Main --------------------------------------------------------------------

cutoff <- 1e-5

pol2.filtered <- LoadPol2Activities()  # $signal and $sitecounts
pol2.signal.all <- LoadPol2Activities.all()
load("/home/yeung/projects/tissue-specificity/results/fits.Robj")  # fit.list from find.oscillating.genes.R to filter rhythmic gene

pol2.filtered$signal <- log2(pol2.filtered$signal + 1)
pol2.signal.all <- log2(pol2.signal.all + 1)

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
dat.liver.rhythmic <- rbind(dat.liver.mrna, dat.liver.pol2)
head(dat.liver.rhythmic)

PlotAmpPhase(dat.liver.rhythmic)


