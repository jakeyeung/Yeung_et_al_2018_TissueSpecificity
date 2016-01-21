# 2015-11-19
# Jake Yeung
# Quantify significance of sitecounts.

library(ggplot2)
library(reshape2)
library(dplyr)
library(hash)

dist.ref <- 500  # 500 left and right of promoter is reference
# sitecounts.path <- "/home/yeung/projects/tissue-specificity/data/sitecounts/motevo/sitecount_matrix_geneids"
# dist <- 500
sitecounts.path <- "/home/yeung/projects/tissue-specificity/data/sitecounts/motevo/encode_50000_dist_sum_multigene/sitecounts.50000.multigene.mat"
dist <- 50000  # needs rescaling
# sitecounts.path <- "/home/yeung/projects/tissue-specificity/data/sitecounts/motevo/encode_1000_dist_sum_multigene/sitecounts.1000.multigene.mat"
# dist <- 1000

# Function ----------------------------------------------------------------


# Load --------------------------------------------------------------------

load("Robjs/dat.long.fixed_rik_genes.Robj")
load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.15.phase_sd_maxdiff_avg.Robj")
load("Robjs/N.long.promoters_500.Robj")

source("scripts/functions/FisherTestSitecounts.R")
source("scripts/functions/PlotGeneAcrossTissues.R")

# Compare with promoters --------------------------------------------------


# models.tw <- sort(as.character(unique(subset(fits.best, n.rhyth >= 8)$model)))
models.flat <- ""

# Tissue-wide
fits.tw <- subset(fits.best, n.rhyth >= 8)
models.tw <- unique(as.character(fits.tw$model))

jmodel <- "Adr"
jtiss <- c("Adr")

jmodel <- "BFAT"
jtiss <- c("BFAT")

jmodel <- "Liver"
jtiss <- c("Liver")

cutoffs <- seq(from = 0.4, to = 0.8, by = 0.1)
N.sub <- RunFisherOnPromoters(N.long, foreground.models = models.tw, background.models = models.flat, cutoffs = cutoffs)



# fits.adrbfataorta <- subset(fits.best, n.rhyth == 3)
# fits.adrbfataorta <- fits.adrbfataorta[grep("Adr.*Aorta.*BFAT", fits.adrbfataorta$model), ]
# length(unique(fits.adrbfataorta$gene))
# jmodel <- unique(as.character(fits.adrbfataorta$model))
# jtiss <- jmodel
# 
# fits.bfataortamus <- subset(fits.best, n.rhyth == 3)
# fits.bfataortamus <- fits.bfataortamus[grep("Aorta.*BFAT.*Mus", fits.bfataortamus$model), ]
# jmodel <- unique(as.character(fits.adrbfataorta$model))
# jtiss <- jmodel
# 
# fits.bfataorta <- subset(fits.best, n.rhyth == 2)
# fits.bfataorta <- fits.bfataorta[grep("Aorta.*BFAT|BFAT.*Aorta", fits.bfataorta$model), ]
# jmodel <- unique(as.character(fits.bfataorta$model))
# # jmodel <- c(jmodel, "BFAT")  # include tissue-specific module BFAT if you want
# jtiss <- jmodel
# 
# fits.livkid <- subset(fits.best, n.rhyth == 2)
# fits.livkid <- fits.livkid[grep("Liver.*Kidney|Kidney.*Liver", fits.livkid$model), ]
# jmodel <- unique(as.character(fits.livkid$model))
# jtiss <- jmodel
# 
# # the Aorta,BFAT antiphasic module
# fits.bfataorta <- subset(fits.best, n.rhyth > 1 & n.rhyth < 11)
# fits.bfataorta <- fits.bfataorta[grep("(;|^)Aorta.*;BFAT(;|$)", fits.bfataorta$model), ]
# jmodel <- unique(as.character(fits.bfataorta$model))
# jtiss <- jmodel
# 
# # Tissue-wide
# fits.tw <- subset(fits.best, n.rhyth >= 8)
# jmodel <- unique(as.character(fits.tw$model))
# jtiss <- jmodel


N.sub.base.livertw <- subset(N.long, model %in% c(models.tw, jmodel))

length(unique(subset(N.sub.base.livertw, model %in% models.tw)$gene))
# length(fits.livkid <- fits.livkid[grep("Liver.*Kidney|Kidney.*Liver", fits.livkid$model), ]$gene)

N.sub.base.livertw$model <- sapply(N.sub.base.livertw$model, function(m){
  if (!m %in% jtiss){
    return("Flat")
  } else {
    return("Rhyth")
  }
})

length(unique(subset(N.sub.base.livertw, model == "Rhyth")$gene))

start <- Sys.time()
cutoffs <- seq(from = 0.4, to = 0.8, by = 0.1)
N.ftest.ltw.all <- data.frame()
for (cutoff in cutoffs){
  print(cutoff)
  N.ltw.ftest <- N.sub.base.livertw %>%
    group_by(motif) %>%
    do(FisherTestSitecounts(., cutoff))
  N.ltw.ftest$cutoff <- cutoff
  N.ftest.ltw.all <- rbind(N.ftest.ltw.all, N.ltw.ftest)
}
print(Sys.time() - start)

N.ftest.ltw.sum <- N.ftest.ltw.all %>%
  group_by(motif) %>%
  summarise(odds.ratio = mean(odds.ratio), p.value = mean(p.value))

ggplot(N.ftest.ltw.sum, aes(y = -log10(p.value), x = odds.ratio, label = motif)) + geom_point() + geom_text()

jcutoff <- 0.6
FisherTestSitecounts(subset(N.sub.base.livertw, motif == "HNF1A.p2"), cutoff=jcutoff, show.table = TRUE)
FisherTestSitecounts(subset(N.sub.base.livertw, motif == "MEF2.A.B.C.D..p2"), cutoff=jcutoff, show.table = TRUE)
FisherTestSitecounts(subset(N.sub.base.livertw, motif == "NFIL3.p2"), cutoff=jcutoff, show.table = TRUE)
FisherTestSitecounts(subset(N.sub.base.livertw, motif == "RORA.p2"), cutoff=jcutoff, show.table = TRUE)
FisherTestSitecounts(subset(N.sub.base.livertw, motif == "ELK1.4_GABP.A.B1..p3"), cutoff=jcutoff, show.table = TRUE)
FisherTestSitecounts(subset(N.sub.base.livertw, motif == "FOX.C1.C2..p2"), cutoff=jcutoff, show.table = TRUE)
FisherTestSitecounts(subset(N.sub.base.livertw, motif == "bHLH_family.p2"), cutoff=jcutoff, show.table = TRUE)

