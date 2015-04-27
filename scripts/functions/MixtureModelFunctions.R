# MixtureModelFunctions.R
# 2015-04-27

ShannonEntropyMixMdl <- function(x, mixmdl){
  p1 <- PredictFromMM(mixmdl, x)
  p2 <- 1 - p1
  entropy <- 0
  for (p in c(p1, p2)){
    entropy <- entropy + (p * log(1 / p))
  }
  return(entropy)
}


GetDensity <- function(mixmdl, x, i){
  # Get density from a single distribution, given by i
  density.i <- mixmdl$lambda[i] * dnorm(x, mean = mixmdl$mu[i], sd = mixmdl$sigma[i])
}

MixtureDensity <- function(mixmdl, x){
  # http://www.r-bloggers.com/fitting-mixture-distributions-with-the-r-package-mixtools/
  # estimates density kernel from mixture model as a function of x
  # p(x) = lambda[1] n(x; mu[1], sigma[1]) + lambda[2] n(x; mu[2], sigma[2])
  p.of.x <- GetDensity(mixmdl, x, 1) + GetDensity(mixmdl, x, 2)
}

PredictFromMM <- function(mixmdl, x){
  # Given mixture model and a value of x,
  # what is probability that the value came from
  # the FIRST gaussian curve.
  # 
  # Assumes mixture model of two gaussians
  density.1 <- GetDensity(mixmdl, x, 1)
  p.of.x <- MixtureDensity(mixmdl, x)
  p1.given.x <- density.1 / p.of.x
  return(p1.given.x)
}