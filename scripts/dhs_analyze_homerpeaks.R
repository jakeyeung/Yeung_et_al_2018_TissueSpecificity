# Jake Yeung
# 2015-03-26
# dhs_analyze_homerpeaks.R
library(ggplot2)

# Functions ---------------------------------------------------------------

GetPeakwidth <- function(col1){
  # Given column 1 row 1 of autocorrelation.txt from homer,
  # get Peakwidth
  peakwidth <- strsplit(col1, split = '\\.\\.')[[1]][[4]] 
  peakwidth <- as.numeric(peakwidth)
  return(peakwidth)
}

GetFragLength <- function(col1){
  # given col1 row1 get fraglength from autocorrelation.txt
  fraglength <- as.numeric(strsplit(col1, split = '\\.\\.')[[1]][[2]])
  return(fraglength)
}

Cnames <- function(){
  # col2 is distance
  # col2 and col3 just + and - strand
  return(c("distance", "pos", "neg"))
}

Cnames.long <- function(){
  return(c("distance", "n.obs", "strand", "peakwidth", "fraglength", "samp"))
}

ReadAutocorrelation <- function(f.path, samp.name){
  # Read Autocorrelation.txt from homer output
  samp.ac <- read.table(f.path, sep='\t', header=TRUE)
  peakwidth <- GetPeakwidth(colnames(samp.ac)[1])
  fraglength <- GetFragLength(colnames(samp.ac)[1])
  colnames(samp.ac) <- Cnames()  # makes calling vectors easier
  # column 1 header contains peakwidth info
  samp.ac.long <- data.frame(distance = rep(samp.ac$distance, 2),
                             n.obs = c(samp.ac$pos, samp.ac$neg),
                             strand = c(rep("+", length(samp.ac$pos)), rep("-", length(samp.ac$neg))),
                             peakwidth = rep(peakwidth, nrow(samp.ac) * 2),
                             fraglength = rep(fraglength, nrow(samp.ac) * 2),
                             samp = rep(samp.name, nrow(samp.ac) * 2))
  cnames <- Cnames.long()
  colnames(samp.ac.long) <- cnames
  return(samp.ac.long)
}

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
  dat.ac <- ReadAutocorrelation(f.path, samp.base)
  dat.merged <- rbind(dat.merged, dat.ac)
}


# Create dat labels for plotting ------------------------------------------


dat.labels <- ddply(dat.merged, .(samp), function(d){
  sprintf("PeakWidth==%s;FragLength==%s", unique(d$peakwidth), unique(d$fraglength))
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
                     samp %in% suspects), 
       aes(x = distance, y = n.obs, group = strand, colour = strand)) +
  geom_line() +
  facet_wrap(~samp) +
  geom_text(data = subset(dat.labels, samp %in% suspects), 
            aes(x = 0, y = 0, label = jlabel), 
            colour = "black", 
            parse = TRUE)
