FitRhythmic <- function(dat, T = 24){
  # Fit rhythmic model to long dat. Expect dat to be per tissue per gene.
  # use lapply and ddply to vectorize this code.
  # dat needs columns: tissue time experiment exprs
  
  
  # Get parameters from complex model: used later
  GetParamsRhythModel <- function(myfit){
    model.params <- coef(myfit)
    return(list(intercept.array = model.params[1],
                intercept.rnaseq = model.params[2],
                a = model.params[3],
                b = model.params[4]))
  }
  
  tissue <- unique(dat$tissue)
  w = 2 * pi / T
  # Expect columns exprs and time
  rhyth.fit <- lm(exprs ~ 0 + experiment + cos(w * time) + sin(w * time), data = dat)
  flat.fit <- lm(exprs ~ 0 + experiment, data = dat)
  compare.fit <- anova(flat.fit, rhyth.fit)
  pval <- compare.fit["Pr(>F)"][[1]][2]
  model.params <- GetParamsRhythModel(rhyth.fit)  # y = experimentarray + experimentrnaseq + a cos(wt) + b sin(wt)
  amp <- sqrt(model.params$a ^ 2 + model.params$b ^ 2)
  phase.rad <- atan2(model.params$b, model.params$a)
  phase.time <- (phase.rad / w) %% T
  dat.out <- data.frame(tissue = tissue, 
                       cos.part = model.params$a, sin.part = model.params$b, 
                       amp = amp, phase = phase.time, pval = pval,
                       int.array = model.params$intercept.array, int.rnaseq = model.params$intercept.rnaseq)
  return(dat.out)
}

FitRhythmicDatLong <- function(dat.long){
  library(parallel)
  dat.long.by_genetiss <- group_by(dat.long, gene, tissue)
  dat.long.by_genetiss.split <- split(dat.long.by_genetiss, dat.long.by_genetiss$tissue)
  print("Finding rhythmic genes (~3 minutes)")
  start <- Sys.time()
  dat.fitrhyth.split <- mclapply(dat.long.by_genetiss.split, function(jdf){
    do(.data = jdf, FitRhythmic(df = .))
  }, mc.cores = 12)
  dat.fitrhyth <- do.call(rbind, dat.fitrhyth.split)
  print(Sys.time() - start)
  return(dat.fitrhyth)
}

