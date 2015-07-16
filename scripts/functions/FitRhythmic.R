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
    rhyth <- jdf %>%
      group_by(gene) %>%
      do(FitRhythmic(dat = .))
  }, mc.cores = 12)
  dat.fitrhyth <- do.call(rbind, dat.fitrhyth.split)
  print(Sys.time() - start)
  return(dat.fitrhyth)
}

FindMostRhythmic <- function(dat, colname="amp", decreasing = TRUE){
  # return first row after sorting by colname
  return(dat[order(dat[[colname]], decreasing = decreasing), ][1, ])
}

GetAmpRelToMax <- function(dat, fits.mostrhythmic){
  tissue <- as.character(dat$tissue[1])
  amp.max <- fits.mostrhythmic[tissue, ]$amp 
  dat$relamp <- dat$amp / amp.max
  return(dat)
}

GetRelamp <- function(fits, max.pval = 1e-3){
  # using FindMostRhythmic and GetAmpRelToMax together
  fits.mostrhythmic <- fits %>%
    group_by(tissue) %>%
    filter(pval <= max.pval) %>%
    #   filter(pval.adj <= max.pvaladj) %>%
    do(FindMostRhythmic(.)) %>%
    data.frame(.)
  rownames(fits.mostrhythmic) <- fits.mostrhythmic$tissue  # indexing
  
  # Get amplitude relative to max amp ---------------------------------------
  fits.relamp <- fits %>%
    group_by(tissue) %>%
    do(GetAmpRelToMax(., fits.mostrhythmic))
  return(fits.relamp)
}

GetAmpRelToGene <- function(dat, base.amp.dic){
  tissue <- as.character(dat$tissue[1])
  amp.base <- base.amp.dic[[tissue]]
  dat$relamp <- dat$amp / amp.base
  return(dat)
}

GetRelampByGene <- function(fits, by.gene = "Arntl"){
  library(hash)
  fits.base <- subset(fits, gene == by.gene)
  base.amp.dic <- hash(as.character(fits.base$tissue), fits.base$amp)
  fits.relamp <- fits %>%
    group_by(tissue) %>%
    do(GetAmpRelToGene(., base.amp.dic))
  return(fits.relamp)
}

IsRhythmicApply <- function(x, pval.min = 1e-4, pval.max = 0.05, relamp.max = 0.1, cutoff = 5.6, pvali = 7, relampi = 10, avg.exprsi = 9){
  # hard-coded shit... not the best but works fast, you should do quick sanity check 
  pval <- as.numeric(x[pvali])
  relamp <- as.numeric(x[relampi])
  avg.exprs <- as.numeric(x[avg.exprsi])
  if (is.nan(pval) | avg.exprs < cutoff){
    return(NA)
  }
  if (pval <= pval.min & relamp >= relamp.max){
    return(TRUE)
  } else if (pval >= pval.max & relamp <= relamp.max){
    return(FALSE)
  } else {
    return("Between")
  }
}
