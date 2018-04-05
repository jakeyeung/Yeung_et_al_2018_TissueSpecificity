# Jake Yeung
# 2015-03-23
# plot_activities_polar_coordinates.R
# Saeed's beautiful activities need to be plotted in polar coordinates.

library(devtools)
# install_github(repo = "naef-lab/f24")  # IF NECESSARY
library(f24.R2.cycling)

A.path <- "/home/yeung/projects/tissue-specificity/results/ridge_mara_rhythmicgenes_liver_activities.txt"
outpath <- "/home/yeung/projects/tissue-specificity/results/ridge_mara_rhythmicgenes_liver_activities.f24.txt" 
A.mat <- read.table(A.path, header = FALSE, sep = '\t', row.names = 1)

t.start <- 18
t.interval <- 2  # every 2 hours
cutoff <- 1e-8
cutoff.amp <- 0.035
n.timepoints <- ncol(A.mat)
t.end <- (n.timepoints - 1) * t.interval + t.start
t.vec <- seq(t.start, t.end, t.interval)

colnames(A.mat) <- t.vec

x <- A.mat[1, ]

# run F24
Y <- data.frame(t(apply(A.mat, 1, f24_R2_cycling, t.vec)))

# get pval adjust
Y$pval.adj <- p.adjust(Y$pval)

# write to table
write.table(Y, outpath, quote = FALSE, sep = "\t", row.names = TRUE, col.names = NA)

Y$gene <- rownames(Y)
Y$is.rhythmic <- sapply(Y$pval.adj, function(x) x < cutoff)
Y$is.highamp <- sapply(Y$amp, function(x) x > cutoff.amp)

Y$label <- rep("", nrow(Y))  # init
Y$label[which(Y$is.rhythmic | Y$is.highamp)] <- Y$gene[which(Y$is.rhythmic | Y$is.highamp)]

# plot 
ggplot(data = Y, aes(x = amp, y = phase, label = label)) + 
  geom_point() + 
  coord_polar(theta = "y") +
  xlab("Amplitude of activity") +
  ylab("Phase of activity") +
  geom_text(aes(x = amp, y = phase), hjust = 0, vjust = 0) +
  ggtitle(paste0("Labeled if adj.pval cutoff<", cutoff, " amp>", cutoff.amp)) +
  scale_y_continuous(limits = c(0, 24), breaks = seq(2, 24, 2))

