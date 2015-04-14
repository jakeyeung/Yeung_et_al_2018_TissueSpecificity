PlotActivitiesWithSE <- function(df){
  jgene <- unique(df$gene)
  ggplot(df, 
         aes(x = time, y = exprs, group = experiment, colour = experiment)) +
    geom_line() +
    geom_errorbar(aes(ymax = exprs + se, ymin = exprs - se)) +
    facet_wrap(~tissue) + 
    xlab("CT") +
    ylab("Activity") + 
    ggtitle(jgene)
}

PlotMeanActivitiesWithSE <- function(df){
  jgene <- unique(df$gene)
  ggplot(df,
         aes(x = tissue, y = exprs, group = experiment, colour = experiment)) + 
    geom_line() + 
    geom_errorbar(aes(ymax = exprs + se, ymin = exprs - se)) +
    ylab("Activity") +
    ggtitle(jgene)
}

GetExprsMean <- function(df){
  # Input: long df filtered for a tissue and a gene.
  # Expect "exprs", "se" and "experiment"
  df.split <- split(df, df$experiment)
  
  # get exprs
  exprs.array.mean <- mean(df.split$array$exprs)
  exprs.rnaseq.mean <- mean(df.split$rnaseq$exprs)
  
  # get se
  var.array.mean <- mean(df.split$array$se^2)
  se.array.mean <- sqrt(var.array.mean)
  var.rnaseq.mean <- mean(df.split$rnaseq$se^2)
  se.rnaseq.mean <- sqrt(var.rnaseq.mean)
  
  df.out <- data.frame(tissue = unique(df$tissue),
                       exprs = c(exprs.array.mean, exprs.rnaseq.mean),
                       se = c(se.array.mean, se.rnaseq.mean),
                       experiment = c("array", "rnaseq"))
}


IsRnaseq <- function(label.samp){
  # Split by period, if splits into two elements, it's RNASeq, otherwise it's array
  split.length <- length(strsplit(label.samp, "[.]")[[1]])
  if (split.length == 2){
    return(TRUE)
  } else if (split.length == 1){
    return(FALSE)
  }
}

GetMergedColnames <- function(cnames.merged){
  # When we load merged table, our colnames have .1 to represent the RNA-Seq.
  # Fix so that it is Sample.Experiment e.g. Adr18.array or Adr18.rnaseq
  cnames.fixed <- sapply(cnames.merged, function(s){
    # Add .array suffix or .rnaseq suffix depending on if it is Rnaseq or Array
    s.label <- strsplit(s, "[.]")[[1]][[1]]
    if (IsRnaseq(s)){
      s.fixed <- paste0(s.label, ".rnaseq")
    } else {
      s.fixed <- paste0(s.label, ".array")
    }
    return(s.fixed)
  })
}

FitRhythmicWeighted <- function(df, T = 24){
  # Input: long df with exprs, time, se, experiment (array or rnaseq) for
  # a given tissue and a given gene. 
  # 
  # Fits a weighted regression model. With weights in SE
  w = 2 * pi / T  # omega
  tissue <- unique(df$tissue)
  
  sigma.sq <- df$se ^ 2
  jweights <- 1 / sigma.sq
  fit.rhyth <- lm(exprs ~ 0 + experiment + sin(w * time) + cos(w * time), data = df, weights = jweights)
  fit.flat <- lm(exprs ~ 0 + experiment, data = df, weights = jweights)
  ftest <- anova(fit.flat, fit.rhyth)
  ftest.pval <- ftest[["Pr(>F)"]][[2]]
  amp <- unname(sqrt(coef(fit.rhyth)["sin(w * time)"] ^ 2 + coef(fit.rhyth)["cos(w * time)"] ^ 2))
  phase <- unname(atan2(coef(fit.rhyth)["sin(w * time)"], coef(fit.rhyth)["cos(w * time)"]))
  if (phase < 0){
    phase <- phase + 2 * pi
  }
  phase <- phase / w
  df.out <- data.frame(tissue = tissue, amp = amp, phase = phase, pval = ftest.pval)
}
