# 2016-04-24
# Jake Yeung
# summarize_nconds_analysis.R

rm(list=ls())

source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/PlotFunctions.R")

# Load --------------------------------------------------------------------

load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.15.phase_sd_maxdiff_avg.Robj", v=T)
load("Robjs/dat.long.fixed_rik_genes.Robj", v=T)
load("Robjs/dat.complex.fixed_rik_genes.Robj", v=T)


# How do we summarize everything? -----------------------------------------

# just how many genes are in each n.rhyth
fits.summary <- fits.best %>%
  group_by(n.rhyth) %>%
  summarise(n.genes = length(gene))

ggplot(subset(fits.summary, n.rhyth > 0), aes(x = n.rhyth, y = n.genes)) + geom_bar(stat = "identity")

# label gene with n.rhyth
dat.complex.tempvargene <- dat.complex %>%
  group_by(gene) %>%
  summarise(total.temp.var = sum(Mod(2 * exprs.transformed) ^ 2))

NrhythToStr <- function(n.rhyth, max.n.rhyth = 8){
  if (n.rhyth < max.n.rhyth){
    return(as.character(n.rhyth))
  } else {
    return(paste0(max.n.rhyth, "+"))
  }
}

n.rhyth.hash <- hash(as.character(fits.best$gene), sapply(fits.best$n.rhyth, function(n) NrhythToStr(n)))
n.rhyth.hash2 <- hash(as.character(fits.best$gene), fits.best$n.rhyth)

dat.complex.tempvargene$n.rhyth <- sapply(as.character(dat.complex.tempvargene$gene), function(g){
  n.rhyth <- n.rhyth.hash[[g]]
  if (!is.null(n.rhyth)){
    return(n.rhyth)
  } else {
    return(NA)
  }
})
dat.complex.tempvargene$n.rhyth2 <- sapply(as.character(dat.complex.tempvargene$gene), function(g){
  n.rhyth <- n.rhyth.hash2[[g]]
  if (!is.null(n.rhyth)){
    return(n.rhyth)
  } else {
    return(NA)
  }
})

# summarize
dat.complex.tempvargene.sum <- dat.complex.tempvargene %>%
  group_by(n.rhyth) %>%
  summarise(total.temp.var = sum(total.temp.var))
dat.complex.tempvargene.sum2 <- dat.complex.tempvargene %>%
  group_by(n.rhyth2) %>%
  summarise(total.temp.var = sum(total.temp.var))

dat.complex.tempvargene.sum <- OrderDecreasing(dat = dat.complex.tempvargene.sum, jfactor = "n.rhyth", jval = "total.temp.var")
ggplot(dat.complex.tempvargene.sum, aes(x = n.rhyth, y = total.temp.var)) + geom_bar(stat = "identity")
ggplot(dat.complex.tempvargene.sum2, aes(x = n.rhyth2, y = total.temp.var)) + geom_bar(stat = "identity")


# Violins -----------------------------------------------------------------

fits.best$n.rhyth.fac <- as.factor(sapply(as.numeric(fits.best$n.rhyth), function(n) NrhythToStr(n)))

# ggplot(subset(fits.best, n.rhyth != 0), aes(x = as.factor(n.rhyth.fac), y = 2 * amp.avg)) + geom_jitter(alpha = 0.2) + xlab("Number of rhythmic tissues") + ylab("Average amplitude (peak to trough)") 
ggplot(subset(fits.best, n.rhyth != 0), aes(x = as.factor(n.rhyth.fac), y = 2 * amp.avg)) + geom_boxplot() + xlab("Number of rhythmic tissues") + ylab("Average amplitude (peak to trough)")  + theme_bw(16)
# ggplot(subset(fits.best, n.rhyth != 0 & amp.avg > 0.1), aes(x = as.factor(n.rhyth.fac), y = 2 * amp.avg)) + geom_boxplot() + xlab("Number of rhythmic tissues") + ylab("Average amplitude (peak to trough)")  + theme_bw(16)


