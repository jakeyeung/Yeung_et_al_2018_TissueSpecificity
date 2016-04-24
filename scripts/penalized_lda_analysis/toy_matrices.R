# toy example for LDA
# two motifs are important, the rest are noise


# Source ------------------------------------------------------------------

source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/PlotFunctions.R")
source("scripts/functions/ShannonEntropy.R")
source("scripts/functions/SitecountsFunctions.R")
source("scripts/functions/LdaFunctions.R")
source("scripts/functions/RemoveP2Name.R")
source("scripts/functions/PlotUCSC.R")  # CoordToBed
source("scripts/functions/GetTopMotifs.R")
source("scripts/functions/HandleMotifNames.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/SvdFunctions.R")
source("scripts/functions/FitRhythmic.R")
source("scripts/functions/NcondsFunctions.R")
source("scripts/functions/Queue.R")
source("scripts/functions/ListFunctions.R")
source("scripts/functions/MemoryManagement.R")
source("scripts/functions/EnvironmentsFunctions.R")


# Main --------------------------------------------------------------------



set.seed(0)

n.motifs <- 190
nsamps <- 100
mean1 <- 0
mean2 <- 1
sd.all <- 0.1
sd.noise <- 1
m1 <- matrix(data = c(rnorm(nsamps, mean = mean2, sd = sd.all), rnorm(nsamps, mean = mean2, sd = sd.all)), nrow = nsamps, ncol = 2)
m2 <- matrix(data = c(rnorm(nsamps, mean = mean1, sd = sd.all), rnorm(nsamps, mean = mean1, sd = sd.all)), nrow = nsamps, ncol = 2)
m3 <- matrix(data = c(rnorm(nsamps, mean = mean2, sd = sd.all), rnorm(nsamps, mean = mean1, sd = sd.all)), nrow = nsamps, ncol = 2)

m.all <- rbind(m1, m2, m3)
labs <- c(rep(1, nrow(m1)), rep(2, nrow(m2)), rep(3, nrow(m3)))
colnames(m.all) <- c("TissueFactor", "ClockFactor")

plot(m.all[, "TissueFactor"], m.all[, "ClockFactor"], pch = ".")
text(m.all[, "TissueFactor"], m.all[, "ClockFactor"], labels = labs)

m.noise <- matrix(NA, nrow = nrow(m.all), ncol = n.motifs - ncol(m.all))

mean.noise <- 0.5
sd.noise <- 1
for (i in seq(ncol(m.noise))){
  m.noise[, i] <- matrix(data = rnorm(n = nrow(m.all), mean = mean.noise, sd = sd.noise))
}

colnames(m.noise) <- paste("NoiseFactor", seq(ncol(m.noise)), sep = "_")

m.all.withnoise <- cbind(m.all, m.noise)

plot(m.all.withnoise[, 4], m.all.withnoise[, 10], col = labs)
plot(m.all[, "TissueFactor"], m.all[, "ClockFactor"], pch = ".")
text(m.all[, "TissueFactor"], m.all[, "ClockFactor"], labels = labs)


# Do pLDA -----------------------------------------------------------------

m.out3 <- PenalizedLDA(m.all.withnoise, labs, lambda = 0.01, K = 2, standardized = FALSE)

PlotLdaOut(m.out3)

print(qplot(m.out3$xproj[, 1], m.out3$xproj[, 2], colour = as.factor(labs), geom = "point", alpha = I(1)) + ggtitle("2D plot single factors"))

PlotLdaOut(m.out3, jdim = 1, jtitle = "Discrim 1: single factors")
PlotLdaOut(m.out3, jdim = 2, jtitle = "Discrim 2: single factors")

PlotLdaOut2D(m.out3, jcex = 0.5)

print(length(m.out3))


# Do cross products -------------------------------------------------------

# do crosses
m.all.withnoise.crosses <- CrossProduct(as.data.frame(m.all.withnoise), remove.duplicates = TRUE)
# add single factors
m.all.withnoise.cross.all <- cbind(m.all.withnoise, m.all.withnoise.crosses)
# remove columns with 0 variance 
m.all.withnoise.cross.all[which(colSums(m.all.withnoise.cross.all) == 0)] <- list(NULL)
dim(m.all.withnoise.cross.all)

convert3to2 <- FALSE
distfilt <- 5000

# jlambda <- 0.029  # kidliv
jlambda <- 0.015  # liv only
if (convert3to2){
  labs.3to2 <- labs; labs.3to2[which(labs.3to2 == 3)] <- 2
  out.cross <- PenalizedLDA(m.all.withnoise.cross.all, labs.3to2, lambda = jlambda, K = 1, standardized = FALSE)
  m <- SortLda(out.cross)
  print(length(m))
  BoxplotLdaOut(out.cross, jtitle = "Cross product separation")
  PlotLdaOut(out.cross, take.n = 50, from.bottom = TRUE, jtitle = paste0("Cross product loadings (from bottom). Dist:", distfilt, "\nN FG peaks:", length(unique(mat.fg$peak)), "\nN BG peaks:", length(unique(mat.bgnonliver$peak))))
  PlotLdaOut(out.cross, take.n = 50, from.bottom = FALSE, jtitle = paste0("Cross product loadings (from top). Dist:", distfilt, "\nN FG peaks:", length(unique(mat.fg$peak)), "\nN BG peaks:", length(unique(mat.bgnonliver$peak))))
} else {
  out.cross <- PenalizedLDA(m.all.withnoise.cross.all, labs, lambda = jlambda, K = 2, standardized = FALSE)
  print(qplot(out.cross$xproj[, 1], out.cross$xproj[, 2], colour = as.factor(labs), geom = "point", alpha = I(1)) + ggtitle("2D plot single factors"))
  plot(out.cross$xproj[, 1], out.cross$xproj[, 2], pch = ".")
  text(out.cross$xproj[, 1], out.cross$xproj[, 2], labels = labs)
  
  indx <- grepl("^TissueFactor$|^ClockFactor$|TissueFactor;ClockFactor", names(out.cross$x))
  factornames <- names(out.cross$x)[indx]
  discrim1 <- out.cross$discrim[indx, 1]; discrim2 <- out.cross$discrim[indx, 2]
  plot(discrim1, discrim2, pch = ".")
  text(discrim1, discrim2, labels = factornames)
  abline(h = 0); abline(v = 0)
}


# Show that a single dimension can separate -------------------------------

plot(m.all.withnoise.cross.all[, "TissueFactor;ClockFactor"], m.all.withnoise.cross.all[, "NoiseFactor_3"], pch = ".")
text(m.all.withnoise.cross.all[, "TissueFactor;ClockFactor"], m.all.withnoise.cross.all[, "NoiseFactor_3"], labels = labs)

plot(m.all.withnoise.cross.all[, "TissueFactor"], m.all.withnoise.cross.all[, "NoiseFactor_3"], pch = ".")
text(m.all.withnoise.cross.all[, "TissueFactor"], m.all.withnoise.cross.all[, "NoiseFactor_3"], labels = labs)

plot(m.all.withnoise.cross.all[, "ClockFactor"], m.all.withnoise.cross.all[, "NoiseFactor_3"], pch = ".")
text(m.all.withnoise.cross.all[, "ClockFactor"], m.all.withnoise.cross.all[, "NoiseFactor_3"], labels = labs)

# Do cross product on rhyth and liver only --------------------------------

mat.clockfactor <- subset(as.data.frame(m.all.withnoise), select = "ClockFactor")
mat.tissfactor <- subset(as.data.frame(m.all.withnoise), select = "TissueFactor")
mat.tissclock <- CrossProductTwoSets(mat.clockfactor, mat.tissfactor)
out.cross <- PenalizedLDA(cbind(m.all.withnoise, mat.tissclock), labs.3to2, lambda = jlambda, K = 1, standardized = FALSE)
m <- SortLda(out.cross)
print(length(m))
BoxplotLdaOut(out.cross, jtitle = "Cross product separation")
PlotLdaOut(out.cross, take.n = 50, from.bottom = TRUE, jtitle = paste0("Cross product loadings (from bottom). Dist:", distfilt, "\nN FG peaks:", length(unique(mat.fg$peak)), "\nN BG peaks:", length(unique(mat.bgnonliver$peak))))
PlotLdaOut(out.cross, take.n = 50, from.bottom = FALSE, jtitle = paste0("Cross product loadings (from top). Dist:", distfilt, "\nN FG peaks:", length(unique(mat.fg$peak)), "\nN BG peaks:", length(unique(mat.bgnonliver$peak))))
