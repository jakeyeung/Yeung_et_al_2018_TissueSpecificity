# Jake Yeung
# 2015-04-17
# Find rhythmic splicing

library(plyr)
library(doParallel)
library(ggplot2)

# Functions ---------------------------------------------------------------
# source("~/projects/tissue-specificity/scripts/miso_analysis/miso_functions/LoadMisoSummary.R")
source("scripts/functions/MakeCluster.R")
source("scripts/miso_analysis/miso_functions//MisoAnalysisFunctions.R")

# Main --------------------------------------------------------------------

# Make Cluster ------------------------------------------------------------

MakeCluster()


# Define stuff ------------------------------------------------------------


path.se <- "data/miso/SE.filelist.10.union"
summary.se <- LoadMisoSummary(path.se, colname_index = 9)
print(colnames(summary.se))

tissues <- GetTissueMiso(colnames(summary.se))
times <- GetTimeMiso(colnames(summary.se))
print(tissues)
print(times)


# Transform to my thing ---------------------------------------------------

eps <- 1e-2
summary.se.transformed <- log2(1 / (1 - (summary.se - eps)))
range(summary.se.transformed)
plot(density(unlist(summary.se.transformed)))

# Make long ---------------------------------------------------------------

summary.long <- MakeLong(df = summary.se)
summary.long.transformed <- MakeLong(df = summary.se.transformed)

# Find rhythmic -----------------------------------------------------------

summary.split <- split(summary.long, summary.long$tissue)

starttime <- Sys.time()
summary.split.fit <- lapply(summary.split, function(df.tiss){
  ddply(df.tiss, .(miso_id), GetRhythmicMiso, .parallel = TRUE)
})
print(Sys.time() - starttime)

head(summary.split.fit$Liv)

summary.fit <- do.call(rbind, summary.split.fit)

# Find significant pvals --------------------------------------------------

summary.fit$pval.adj <- p.adjust(summary.fit$pval)
head(summary.fit[order(summary.fit$pval.adj), ], n = 20)
head(summary.fit[order(summary.fit$amp, decreasing = TRUE), ], n = 20)

cutoff.pval <- 0.01
cutoff.amp <- 0.1

summary.fit.filter <- subset(summary.fit, pval <= cutoff.pval & amp >= cutoff.amp)
head(summary.fit.filter)

# plot significant miso_ids
jid <- "chr11:115419892:115419962:-@chr11:115419192:115419293:-@chr11:115418386:115418516:-"
jid <- "chr12:81532716:81532907:-@chr12:81510829:81510965:-@chr12:81504544:81504639:-"
jid <- "chr2:121457009:121457291:+@chr2:121457646:121457701:+@chr2:121457992:121458043:+"
jid <- "chr18:35598667:35598759:+@chr18:35609593:35609673:+@chr18:35610871:35611034:+"
jid <- "chr11:75679259:75679664:+@chr11:75692197:75692732:+@chr11:75703366:75708428:+"
jid <- "chr6:29366917:29366977:+@chr6:29372471:29372670:+@chr6:29374399:29376675:+"

jid <- "chr14:52103889:52104028:-@chr14:52098015:52098040:-@chr14:52084115:52084391:-"

PlotMisoAcrossTissue(subset(summary.long, miso_id == jid))

