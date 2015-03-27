# Jake Yeung
# 2015-03-26
# dhs_analyze_homerpeaks.R
library(ggplot2)

# Functions ---------------------------------------------------------------

source("scripts/functions/DhsFunctions.R")

# Main --------------------------------------------------------------------

dat.dir <- "/home/yeung/projects/tissue-specificity/data/homerpeaks"
dir.vec <- list.files(dat.dir)
samp.paths <- sapply(dir.vec, function(d){
  file.path(dat.dir, d)
}, USE.NAMES=FALSE)

# init merged dataframe and peakwidths
cnames <- Cnames.long()
dat.merged <- data.frame(matrix(nrow = 0, ncol = length(cnames)))
colnames(dat.merged) <- cnames

# Loop through files, create merged data in long format -------------------


fname <- "tagAutocorrelation.txt"
for (samp in samp.paths){
  samp.base <- basename(samp)
  f.path <- file.path(samp, fname)
  dat.temp <- ReadAutocorrelation(f.path, samp.base)
  dat.merged <- rbind(dat.merged, dat.temp)
}

# Create dat labels for plotting ------------------------------------------


dat.labels <- ddply(dat.merged, .(samp), function(d){
  sprintf("FragLength==%s;PeakWidth==%s", unique(d$peakwidth), unique(d$fraglength))
  # sprintf("FragLength==%s,Peakwidth==%s", unique(d$peakwidth), unique(d$fraglength))
})
colnames(dat.labels) <- c("samp", "jlabel")
dat.labels$strand = rep("+", nrow(dat.labels))  # colour is strand, messes up the code otherwise


# Plot all samples --------------------------------------------------------

ggplot(data = dat.merged, 
       aes(x = distance, y = n.obs, group = strand, colour = strand)) +
  geom_line() +
  facet_wrap(~samp) + 
  geom_text(data = dat.labels,
            aes(x = 0, y = 0, label = jlabel), 
            colour = "black", 
            parse = TRUE) + 
  scale_x_continuous(limits = c(-200, 200))


# Plot best one -----------------------------------------------------------

best <- c("wgEncodeUwDnaseLiverC57bl6MAdult8wksRawDataRep9")
ggplot(data = subset(dat.merged, 
                     samp %in% best), 
       aes(x = distance, y = n.obs, group = strand, colour = strand)) +
  geom_line() +
  facet_wrap(~samp) +
  geom_text(data = subset(dat.labels, samp %in% best), 
            aes(x = 0, y = 0, label = jlabel), 
            colour = "black", 
            parse = TRUE)

# Plot the smallest two ---------------------------------------------------


suspects <- c("UwStam_mLiver-DS16858-FC62FJY-4_004", 
              "UwStam_mLiver-DS19375-FCC05B7-L006_R1_002")

ggplot(data = subset(dat.merged, 
                     samp %in% suspects), 
       aes(x = distance, y = n.obs, group = strand, colour = strand)) +
  geom_line() +
  facet_wrap(~samp) +
  geom_text(data = subset(dat.labels, samp %in% suspects), 
            aes(x = 0, y = 0, label = jlabel), 
            colour = "black", 
            parse = TRUE)


# Plot the 3 semi-outliers ------------------------------------------------


suspects2 <- c("UwStam_mLiver-DS16853-FC62FJY-3_004",
               "UwStam-mLiver-DS16740-FC62J2C-7-003",
               "UwStam_mLiver-DS19636-FCD06E9-L006_R1_002")

ggplot(data = subset(dat.merged, 
                     samp %in% suspects2), 
       aes(x = distance, y = n.obs, group = strand, colour = strand)) +
  geom_line() +
  facet_wrap(~samp) +
  geom_text(data = subset(dat.labels, samp %in% suspects2), 
            aes(x = 0, y = 0, label = jlabel), 
            colour = "black", 
            parse = TRUE)


# Now do the same for tag distributions -----------------------------------

tcd.cnames <- GetTcdCnames() 
dat.tcd <- data.frame(matrix(nrow = 0, ncol = length(tcd.cnames)))
colnames(dat.tcd)

countdist.fname <- "tagCountDistribution.txt"
for (samp in samp.paths){
  countdist.base <- basename(samp)
  countdist.path <- file.path(samp, countdist.fname)
  dat.temp <- ReadTagCountsDistribution(countdist.path, countdist.base)
  dat.tcd <- rbind(dat.tcd, dat.temp)
}

ggplot(data = dat.tcd, aes(x = tags.per.pos, y = frac.pos)) +
  geom_bar(stat = "identity") +
  scale_x_continuous(breaks = 1:10, limits = c(0, 10)) +
  facet_wrap(~samp)
