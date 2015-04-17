# Jake Yeung
# 2015-04-17
# Find rhythmic splicing

library(plyr)
library(doParallel)
library(ggplot2)

# Functions ---------------------------------------------------------------
# source("~/projects/tissue-specificity/scripts/miso_analysis/miso_functions/LoadMisoSummary.R")
source("scripts/functions/MakeCluster.R")

LoadMisoSummary <- function(path, colname_index=9){
  # Load miso summary file. Fix the column names which are filepaths
  df <- read.table(path)
  if (colname_index != FALSE){
    new_colnames <- sapply(colnames(df), function(x){
      # Take the colname_index from original colname which looks like:
      # X.scratch.el.monthly.jyeung.hogenesch_rna_seq.miso_run_outputs_v2.SE.Adr_CT22Aligned.sortedByCoord.out.bam.summary.summary.Adr_CT22Aligned.sortedByCoord.out.bam.miso_summary.10.filtered
      # Example: colname_index = 9 extracts Adr_CT22Aligned
      fixed_colname <- strsplit(x, "\\.")[[1]][[colname_index]]
      # Get Adr_CT22 from Adr_CT22Aligned
      fixed_colname <- substr(fixed_colname, 1, nchar(fixed_colname) - 7)
      return(fixed_colname)
    })
    colnames(df) <- new_colnames
  }
  return(df)
}

GetTissueMiso <- function(cnames){
  # input: vector f colnames like Adr_CT22
  # Extract Adr from Adr_CT22, return as vector
  tissues <- rep(NA, length(cnames))
  for (i in 1:length(cnames)){
    cname <- cnames[i]
    tissue <- strsplit(cname, '_')[[1]][[1]]
    tissues[i] <- tissue
  }
  return(tissues)
}

GetTimeMiso <- function(cnames, time_prefix="CT"){
  # input vector of colnames like Adr_CT22
  # extract 22 return as vector
  times <- rep(NA, length(cnames))
  for (i in 1:length(cnames)){
    cname <- cnames[i]
    time <- strsplit(cname, '_')[[1]][[2]]
    # remove CT from CT22
    time <- strsplit(time, time_prefix)[[1]][[2]]
    times[i] <- time
}
  return(as.numeric(times))
}

GetRhythmicMiso <- function(df, T = 24){
  # Expect miso_id, psi, tissue, time column names
  # df is for one tissue and one miso_id
  
  # Define omega
  w <- 2 * pi / T
  
  # Fit rhythmic and flat, then test for rhythmicity
  fit.rhyth <- lm(psi ~ sin(w * time) + cos(w * time), data = df)
  fit.flat <- lm(psi ~ 1, data = df)
  test.rhyth <- anova(fit.flat, fit.rhyth)
  test.pval <- test.rhyth[["Pr(>F)"]][[2]]
  
  # Define amp and phase
  sin.part  <- fit.rhyth$coefficients[["sin(w * time)"]]
  cos.part <- fit.rhyth$coefficients[["cos(w * time)"]]
  intercept.part <- fit.rhyth$coefficients[["(Intercept)"]]
  amp <- sqrt(sin.part ^ 2 + cos.part ^ 2)
  phase <- atan2(sin.part, cos.part)
  
  # output in format that will work well for ddply
  df.out <- data.frame(miso_id = unique(df$miso_id),
                       tissue = unique(df$tissue),
                       amp = amp,
                       phase = phase,
                       intercept = intercept.part,
                       pval = test.pval)
}

PlotMisoAcrossTissue <- function(df){
  # take subset df of a miso_id across tissues and time,
  # plot rhythmicity over time
  # 
  # Expect psi, tissue, time and miso_id as colnames
  m <- ggplot(data = df, aes(x = time, y = psi)) +
              geom_point() +
              geom_line() +
              facet_wrap(~tissue)
  return(m)
}

MakeLong <- function(df){
  df.long <- data.frame(miso_id = rep(rownames(df), ncol(df)),
                        psi = unlist(df),
                        tissue = rep(tissues, each = nrow(df)),
                        time = rep(times, each = nrow(df)))
  return(df.long)
}

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

