# 2015-09-07
# Jake Yeung
# Check that merged data looks OK after cutoff and concatenation

matfile <- "/home/yeung/projects/tissue-specificity/data/beds/merge/cutoffs_smoother_gamma_cutoff_high_accuracy_ALL_gamma_inits_scale50_mu3_12_zeros_normalized_scaled_billion/merged_counts.mat"

dat <- read.table(matfile, header = TRUE, sep = "\t")

dat.bedinfo <- dat[, c("chromo", "start", "end")]

dat.counts <- dat[, 4:ncol(dat)]

hc <- hclust(dist(t(dat.counts)))

plot(hc)

dat.counts.long <- data.frame(tissue = rep(colnames(dat.counts), each = nrow(dat.counts)),
                              count = unlist(dat.counts))

ggplot(dat.counts.long, aes(x = log2(count))) + geom_density() + facet_wrap(~tissue)
