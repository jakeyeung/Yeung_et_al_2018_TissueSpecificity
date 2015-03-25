# Get parameters from complex model
GetParamsRhythModel <- function(myfit){
  model.params <- coef(myfit)
  return(list(intercept.array = model.params[1],
              intercept.rnaseq = model.params[2],
              a = model.params[3],
              b = model.params[4]))
}

FitRhythmic <- function(df, T = 24){
  w = 2 * pi / T
  # Expect columns exprs and time
  rhyth.fit <- lm(exprs ~ 0 + experiment + cos(w * time) + sin(w * time), data = df)
  flat.fit <- lm(exprs ~ 0 + experiment, data = df)
  compare.fit <- anova(flat.fit, rhyth.fit)
  pval <- compare.fit["Pr(>F)"][[1]][2]
  model.params <- GetParamsRhythModel(rhyth.fit)  # y = experimentarray + experimentrnaseq + a cos(wt) + b sin(wt)
  amp <- sqrt(model.params$a ^ 2 + model.params$b ^ 2)
  phase.rad <- atan2(model.params$b, model.params$a)
  phase.time <- (phase.rad / w) %% T
  df.out <- data.frame(cos.part = model.params$a, sin.part = model.params$b, 
                       amp = amp, phase = phase.time, pval = pval,
                       int.array = model.params$intercept.array, int.rnaseq = model.params$intercept.rnaseq)
  return(df.out)
}