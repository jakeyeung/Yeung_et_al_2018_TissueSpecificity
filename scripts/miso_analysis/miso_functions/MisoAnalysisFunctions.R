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
