# find_alternative_first_exon.R
# Analyze MISO alternative first exons across tissues.

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

path.afe <- "data/miso/AFE.mm10.filelist.10.union"
outfile <- "data/miso/AFE.mm10.sighits"

summary.afe <- LoadMisoSummary(path.afe, colname_index = 10)

tissues <- GetTissueMiso(colnames(summary.afe))
times <- GetTimeMiso(colnames(summary.afe))
print(tissues)
print(times)


# Make long ---------------------------------------------------------------

summary.long <- MakeLong(df = summary.afe)


# Mean across tissue ------------------------------------------------------

summary.mean <- ddply(summary.long, .(tissue, miso_id), summarise,
                      Mean = mean(psi),
                      SD = sd(psi))

summary.mean.mean <- ddply(summary.mean, .(miso_id), summarise,
                           Range = diff(range(Mean)))


# Show top between Kidney and Liver ---------------------------------------

sig.hits <- subset(summary.mean.mean, Range >= 0.5)

sink(outfile)
for (id in sig.hits$miso_id){
  cat(id)
  cat("\n")
}
sink()
