# Jake Yeung
# Date of Creation: 2018-04-20
# File: ~/projects/tissue-specificity/scripts/nconds_bayesfactors/milou_test.R
# Milou test C. elegans data

rm(list=ls())

library(ggplot2)
library(f24.R2.cycling)
library(devtools)

install_local("/home/yeung/projects/nconds2")
library(nconds2)  # USB key of nconds_analysis 

source("/home/yeung/projects/nconds_analysis/scripts/functions/ComplexSvdFunctions.R")
# source("/home/yeung/projects/nconds_analysis/scripts/functions/MaraFunctions/MaraDownstream.R")
source("/home/yeung/projects/tissue-specificity/scripts/functions/PlotFunctions.R")

OrderDecreasing <- function(dat, jfactor, jval){
  # Reorder factors by decreasing value of jval
  # used in conjunction with barplots makes life easier
  dat[[jfactor]] <- factor(dat[[jfactor]], levels = dat[[jfactor]][order(dat[[jval]], decreasing = TRUE)])
  return(dat)
}

PlotGene <- function(dat.sub){
  m <- ggplot(dat.sub, aes(x = time, y = exprs, group = tissue, color = tissue)) + 
    geom_point() + geom_line()  + 
    theme_bw()
  return(m)
}

# Load  -------------------------------------------------------------------

inf <- "/home/shared/Milou/mRNAexpression.tab"
inf.chip <- "/home/shared/Milou/polII_TSS_1kb_quantile_norm.tab"

dat <- read.table(inf)
dat <- dat[, -1]
dat.norm <- sweep(dat, MARGIN = 2, STATS = colSums(dat), FUN = "/") * 10^6

dat.chipseq <- read.table(inf.chip)

eps <- 1
dat.log2 <- log2(dat.norm + eps)
dat.means <- rowMeans(dat.log2)
low.exprs.cutoff <- 2.5
lowly.exprs.genes <- names(dat.means[dat.means < low.exprs.cutoff])

jgene <- "WBGene00000599"
plot(seq(ncol(dat.norm)), dat.norm[jgene, ], 'o')

par(mfrow=c(2, 1))
plot(seq(ncol(dat.log2)), dat.log2[jgene, ], 'o')
plot(seq(ncol(dat.chipseq)), dat.chipseq[jgene, ], 'o')
par(mfrow(c(1, 1)))

tvec <- seq(ncol(dat.log2))

fits8 <- data.frame(t(apply(dat.log2, 1, function(row) f24_R2_cycling(x = row, t = tvec, period = 8))))
fits8.chipseq <- data.frame(t(apply(dat.log2, 1, function(row) f24_R2_cycling(x = row, t = tvec, period = 8))))

# Make long ---------------------------------------------------------------

dat.long <- data.frame(gene = rownames(dat.norm),
                       time = rep(tvec, each = nrow(dat.norm)),
                       exprs = unlist(dat.log2),
                       tissue = "mRNA")
dat.long.chipseq <- data.frame(gene = rownames(dat.chipseq),
                               time = rep(tvec, each = nrow(dat.chipseq)),
                               exprs = unlist(dat.chipseq),
                               tissue = "ChIPseq")
dat.merged <- rbind(dat.long, dat.long.chipseq)

ggplot(subset(dat.merged, gene == jgene), aes(x = time, y = exprs, group = tissue, color = tissue)) + 
  geom_point() + geom_line()  + 
  theme_bw() 

dat.env <- DatLongToEnvironment(dat.merged)
tissues.uniq <- as.character(unique(dat.merged$tissue))
tperiod <- 8
test <- MakeDesMatRunFitEnv(dat.env, jgene, tissues.uniq, 
                            n.rhyth.max = 2, w = 2 * pi / tperiod, criterion = "BIC", 
                            normalize.weights = TRUE, cutoff = 1e-5, top.n = NULL, sparse = FALSE)

# Show fits
par(mfrow(c(1, 1)), pty = "s")
plot(seq(ncol(dat.log2)), dat.log2[jgene, ], 'o', col = 'blue')
lines(seq(ncol(dat.chipseq)), dat.chipseq[jgene, ], 'o', col = 'red')

ggplot(subset(dat.merged, gene == jgene), aes(x = time, y = exprs, group = tissue, color = tissue)) + 
  geom_point() + geom_line()  + 
  theme_bw()


# do for real: about 30 minutes
start <- Sys.time()
outf <- "/home/shared/Milou/fit_outputs.Rdata"

print(paste("Outf:", outf))
fits.all <- lapply(ls(dat.env), function(gene){
  MakeDesMatRunFitEnv(dat.env, gene, tissues.uniq,
                      n.rhyth.max = 2, w = 2 * pi / tperiod,
                      criterion = "BIC", normalize.weights = TRUE,
                      cutoff = 1e-5, top.n = NULL, sparse = FALSE)
})
print(Sys.time() - start)

fits.all.long <- lapply(fits.all, function(x){
  gene <- x$gene
  x$gene <- NULL
  fits.temp.long <- ListToLong(x, gene, top.n = 5, period = tperiod)
})
fits.all.long <- do.call(rbind, fits.all.long)
save(fits.all, fits.all.long, file=outf)

print(Sys.time() - start)


# Get best model ----------------------------------------------------------

fits.all.best <- fits.all.long %>%
  group_by(gene) %>%
  filter(weight == max(weight)) %>%
  arrange(desc(weight))

fits.summary <- fits.all.best %>%
  filter(! gene %in% lowly.exprs.genes) %>%
  group_by(model) %>%
  summarise(counts = length(gene))

fits.summary <- OrderDecreasing(fits.summary, "model", "counts")

m.summary <- ggplot(fits.summary, aes(x = model, y = counts)) + geom_bar(stat = "identity") + theme_bw() 

# Find number of clusters -------------------------------------------------

print(head(subset(fits.all.best, model == "ChIPseq")))
jgene <- "WBGene00202164"
PlotGene(subset(dat.merged, gene == jgene))
subset(fits.all.best, model == "ChIPseq")$param.list[[1]]

print(head(subset(fits.all.best, model == "mRNA")))
jgene <- "WBGene00019072"
PlotGene(subset(dat.merged, gene == jgene))
subset(fits.all.best, model == "mRNA")$param.list[[1]]

print(head(subset(fits.all.best, model == "mRNA,ChIPseq")))
jgene <- "WBGene00010045"
PlotGene(subset(dat.merged, gene == jgene))
subset(fits.all.best, model == "mRNA,ChIPseq")$param.list[[1]]

print(head(subset(fits.all.best, model == "mRNA;ChIPseq")))
jgene <- "WBGene00007142"
PlotGene(subset(dat.merged, gene == jgene))
subset(fits.all.best, model == "mRNA;ChIPseq")$param.list[[1]]


plot(density(unlist(dat.chipseq)))



# SVD to summarize --------------------------------------------------------

library(reshape2)
outrobjf <- "/home/shared/Milou/dat_complex.Rdata"
if (!file.exists(outrobjf)){
  # ~ 2 minutes
  dat.complex <- dat.merged %>%
    group_by(gene, tissue) %>%
    do(ProjectToFrequency2(., omega = 2 * pi / tperiod, uneven.samps = TRUE))
  save(dat.complex, file = outrobjf)
}

fits.all.best.filt <- subset(fits.all.best, ! gene %in% lowly.exprs.genes)
  
models.lst <- list(c("mRNA;ChIPseq", "mRNA,ChIPseq"), c("mRNA"), c("ChIPseq"))
pdf("/home/shared/Milou/svd_plot_of_model_selection.pdf")
for (jmod in models.lst){
  genes.sub <- subset(fits.all.best.filt, model %in% jmod)$gene
  
  M.complex <- dcast(subset(dat.complex, gene %in% genes.sub), gene ~ tissue, value.var = "exprs.transformed")
  rownames(M.complex) <- make.names(M.complex$gene, unique = TRUE)
  M.complex$gene <- NULL
  
  s.tw <- SvdOnComplex(M.complex)
  
  comp <- 1
  dotsize <- 3
  
  eigens.tw <- GetEigens(s.tw, period = 8, jsize = 8, comp = comp, label.n = 30, eigenval = TRUE, adj.mag = TRUE, constant.amp = dotsize, peak.to.trough = TRUE)
  multiplot(eigens.tw$u.plot, eigens.tw$v.plot, cols = 2)
  print(eigens.tw$u.plot)
  print(eigens.tw$v.plot)
}
dev.off()


# plot summary of each model
pdf("/home/shared/Milou/model_summary.pdf")
  print(m.summary)
dev.off()

# filtered genes
writeLines(lowly.exprs.genes, con = "/home/shared/Milou/lowly_exprs_genes_removed.txt")

# 
# # Compare ChipSeq with mRNA -----------------------------------------------
# 
# tperiod <- 8
# omega <- 2 * pi / tperiod
# tvec <- seq(1, tperiod + tperiod / 2, length.out = 12)
# # tvec <- seq(1, 24, by = 1)
# jphase <- 0 * omega  # Bmal1
# yvec <- 10 + 2 * cos(omega * tvec - jphase)
# 
# jgene <- "WBGene00007142"
# dat.test <- subset(dat.merged, gene == jgene & tissue == "ChIPseq")
# dat.test <- data.frame(gene = jgene, time = tvec, exprs = yvec)
# plot(dat.test$time, dat.test$exprs, 'o')
# subset(fits.all.best, model == "mRNA;ChIPseq" & gene == jgene)$param.list[[1]]
# 
# (exprs.transformed <- DoFourier(dat.test$exprs, dat.test$time, omega, normalize=FALSE))
# 
# # normalize properly
# tvec.counts <- table(tvec %% tperiod)
# weights <- 1 / tvec.counts[as.character(tvec %% tperiod)]
# weights.norm <- weights / sum(weights)
# exprs.trans <- (weights.norm * dat.test$exprs) %*% exp(1i * omega * dat.test$time)
# 
# # Fourier manually
# dat.test$exprs %*% exp(1i * omega * tvec)
# 
# ConvertArgToPhase(Arg(exprs.transformed), omega)
# GetPhi(a = Re(exprs.transformed), b = Im(exprs.transformed), omega)
# 
# 
# 
# # Check with known data ---------------------------------------------------
# 
# inf <- "/home/yeung/projects/nconds_analysis/Robjs/dat.complex.8conditions.Robj"
# load(inf, v=T)
# 
# jtest <- subset(dat.complex, condition == "RF_Bmal1_WT" & gene == "Arntl")
# ConvertArgToPhase(Arg(jtest$exprs.transformed), 2 * pi / 24)
