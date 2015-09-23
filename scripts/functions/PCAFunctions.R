# used in pca_adjusted_microarray.label_variance.R

PlotComponents <- function(pca.p.med, start, end){
  # replacte gray with some dark blue
  # http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/#a-colorblind-friendly-palette
  cbPalette <- c("#CC79A7", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  break.i <- floor((end - start ) / 4)
  # http://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
  # getting the right colours in different subsets is not trivial!
  m <- ggplot(subset(pca.p.med, pc.num <= end & pc.num >= start), aes(x = pc.num, y = eigenvals, colour = T.max.med, alpha = N / 12)) + 
    geom_bar(stat = "identity", size = 1.0) + scale_x_continuous(breaks=seq(start, end, break.i)) +
    xlab("Component") + ylab("Fraction of total sqr eigenvals") +
    scale_colour_manual(drop=TRUE, limits = levels(pca.p.med$T.max.med), values=cbPalette) +
    scale_alpha_continuous(name = "Frac of tissues \\ with T.max") +
    scale_color_discrete(name = "Mode of T.max")
#   scale_colour_discrete(drop=TRUE, limits = levels(pca.p.med$T.max.med)) 
  print(m)
}

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

PlotPCTissue <- function(dat){
  ggplot(dat, aes(x = time, y = loading)) + geom_point() + geom_line() + facet_wrap(~tissue)
}

SummarisePeriodogram <- function(dat){
  # Summarize periodogram from GetPeriodogramFreq to 
  # get mode T.max, number of tissues with T.max, mean p.pmax for T.max = median(T.max)
  T.max.med <- as.character(Mode(dat$T.max))
  dat.sub <- subset(dat, T.max == T.max.med)
  N <- nrow(dat.sub)
  outdat <- data.frame(T.max.med = T.max.med, N = N)
  return(outdat)
}

GetPeriodogramFreq <- function(dat, interval = 2){
  # interval: hour between each sample. Used to get proper period from frequency
  # 2015-09-14
  x <- dat$loading
  p <- CalculatePeriodogram(x, is.matrix = FALSE)
  max.freq <- FindMaxFreqs(freq = p$freq, periodogram = p$p.scaled, n = 1)
  p.pmax <- p$p.scaled[which(p$freq == max.freq)]
  p.24 <- p$p.scaled[which(p$freq == 2/24)]
  p.inf <- p$p.scaled[which(p$freq == 0)]
  p.time <- sum(p$p.scaled[which(p$freq > 0)])
  # return as dataframe so we can use dplyr
  #   outdat <- list(p.24, p.inf, p.pmax, max.freq)
  outdat <- data.frame(p.24 = p.24, p.inf = p.inf, p.time = p.time, p.pmax = p.pmax, T.max = interval / max.freq)
  return(outdat)
}

MeanCenterAcrossGroups <- function(x, n.per.group=24) {
  # Input:
  # x = vector x of inputs across all samples, contains n.per.group
  # elements per group. 
  # Vector is ordered such that all elements
  # in group i are together. 
  # Counting n.per.group would
  # reveal each group. 
  # 
  # Output:
  # x.mean.center <- calculate mean of each group (i.e. mean
  # of every group of n.per.group) and subtract each element
  # in that group by its mean.
  
  N <- length(x)
  x.mean.center <- rep(NA, length(x))
  # BEGIN SANITY CHECK
  # check that N / n.per.group gives no remainder
  if (N %% n.per.group != 0) {
    warning("Uneven number of elements per group.")
  }
  n.groups <- N / n.per.group  # should be an integer
  # END SANITY CHECK
  
  
  for (i in 0:n.groups - 1) {
    # loop from 0 to n.groups - 1 so that our start index
    # loops like 1, 1 + n.per.group, 1 + 2 * n.per.group
    # for example start.index = 1, 25, 49 ...
    start.index <- i * n.per.group + 1
    end.index <- start.index + n.per.group - 1
    # get group e.g. x[1:24]
    group <- x[start.index:end.index]
    # subtract each element in group by mean(group)
    group.mean.center <- group - mean(group)
    # put into x.mean.center
    x.mean.center[start.index:end.index] <- group.mean.center
  }
  
  return(x.mean.center)
}