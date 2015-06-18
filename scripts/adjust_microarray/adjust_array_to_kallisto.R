# 2015-06-17
# Jake Yeung
# adjust_array_to_kallisto.R
# Adjust array to kallisto gene counts because I trust it more than DESeq

setwd("/home/yeung/projects/tissue-specificity")
library(ggplot2)
library(dplyr)

# Functions ---------------------------------------------------------------

source("scripts/functions/LoadArrayRnaSeq.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/PlotGeneTpm.R")
source("scripts/functions/GetTissueTimes.R")
source("scripts/functions/RegressionFunctions.R")

saturation <- function(params, r) params[2] + (params[1] * r) / (params[3] + r)
saturation.inv <- function(params, m) params[3] * (m - params[2]) / (params[2] + params[1] - m)

saturation2 <- function(a, b, k, r){
  b + (a * r ) / (k + r)
}
saturation.inv2 <- function(a, b, k, m) {
  k * (m - b) / (b + a - m)
}
linear.inv2 <- function(int, slope, y) {
  (y - int) / slope
}


PlotDiagnostics2 <- function(gene, kallisto.wide, array.wide, fits){
  common.samps <- intersect(colnames(kallisto.wide), colnames(array.wide))
  unique.samps <- setdiff(colnames(array.wide), colnames(kallisto.wide))
  
  R <- as.matrix(kallisto.wide[gene, common.samps])
  M <- as.matrix(array.wide[gene, common.samps])
  M.full <- as.matrix(array.wide[gene, ])
  M.uniq <- as.matrix(array.wide[gene, unique.samps])
  
  fit <- fits[gene, ]
  
  if (fit$model == "saturation"){
    yspace <- seq(min(M.full) * 0.9, max(M.full) * 1.1, length.out = 55)
    xspace <- sapply(yspace, function(m) saturation.inv2(a = fit$a, b = fit$b, k = fit$k, m = m))
    
    plot(x = R, y = M, main = gene, xlim = c(min(R), max(R)), ylim = c(min(M.full), max(M.full)))
    lines(xspace, yspace, type = 'l')
    points(x = rep(min(R), length(M.full)), y = M.full, pch = 8, cex = 0.5)
  } else {
    yspace <- seq(min(M.full) * 0.9, max(M.full) * 1.1, length.out = 55)
    xspace <- sapply(yspace, function(m) linear.inv2(int = fit$int, slope = fit$slope, y = m))
    plot(x = R, y = M, main = gene)
    lines(xspace, yspace, type = 'l')
    points(x = rep(min(R), length(M.full)), y = M.full, pch = 8, cex = 0.5)
  }
}

GetGeneList <- function(){
  # for testing
  # Clock genes --------------------------------------------------------
  
  clockgenes <- c('Nr1d1','Dbp', 'Arntl', 'Npas2', 'Nr1d2', 
                  'Bhlhe41', 'Nfil3', 'Cdkn1a', 'Lonrf3', 
                  'Tef', 'Usp2', 'Wee1', 'Dtx4', 'Asb12')
  clockgenes <- c(clockgenes, 'Elovl3', 'Clock', 'Per1', 'Per2', 'Per3', 'Cry2', 'Cry1')
  clockgenes <- c(clockgenes, 'Defb48', 'Svs1', 'Svs2', 'Svs5', 'Defb20', 'Adam7', 'Lcn8', 'Rnase10', 'Teddm1')
  clockgenes <- c(clockgenes, "Sry")  # known to not be expressed via RNASeq
  
  
  # Tissue genes ------------------------------------------------------------
  
  clockgenes <- c(clockgenes, "Plbd1","Acacb","Hadh","Decr1","Acadl","Ech1","Acsl1","Cpt2","Acaa2")  # PCA 1 from explore_data.R
  clockgenes <- c(clockgenes, "Elovl1","Slc35b3","Fkbp15","Slc39a8","Sep15","Mospd2","Med11","Slco2a1","Arf6","Yipf6")  # PCA 2 from explore_data.R
  
  
  # Problem genes -----------------------------------------------------------
  
  # interesting cases
  clockgenes <- c(clockgenes, "Tmem59l", "Saysd1", "H2.Q1", "Col25a1", "Clec18a", "4930427A07Rik", 
                        "Ncan", "Crtam", "Fam43b", "Nphp4", "Nuf2", "Sox11", "Krt23", "Myo1h", 
                        "Mas1", "Cd207", "Tgif2", "Sdsl", "Gm8659", "Fsd1", "X2510049J12Rik", 
                        "1600029I14Rik", "Syndig1l", "Cyp4x1", "E130309D14Rik", "Gm281", 
                        "Wdr27", "Daw1", "Tnfsf9", "Myo16")
  
  # genes with 0 expression
  clockgenes <- c(clockgenes, "1700025F22Rik")
  
  return(clockgenes)
}

GetMeanVarByTissues <- function(exprs, tissue.names){
  # Calculate mean and variance for each gene per tissue
  N <- nrow(exprs) * length(tissue.names)  # one measurement for each gene for all tissues.
  mean.var <- list(mean=matrix(NA, nrow=nrow(exprs), ncol=length(tissue.names),
                               dimnames = list(rownames(exprs), 
                                               tissue.names)), 
                   var=matrix(NA, nrow=nrow(exprs), ncol=length(tissue.names),
                              dimnames = list(rownames(exprs), tissue.names)))
  for (j in 1:length(tissue.names)){
    tissue <- tissue.names[j]
    gene.tissue.exprs <- exprs[, grepl(tissue, colnames(exprs))]
    # calculate mean and var, by row
    exprs.mean <- apply(gene.tissue.exprs, 1, mean)
    exprs.var <- apply(gene.tissue.exprs, 1, var)
    # append to matrix 
    mean.var$mean[, j] <- exprs.mean
    mean.var$var[, j] <- exprs.var
  }
  # Make dataframe for ggplot2 and other analyses
  mean.var.df <- data.frame(mean=as.vector(mean.var$mean), var=as.vector(mean.var$var))
  return(mean.var.df)
}

PlotArrayRnaseq <- function(dat, log2.transform=FALSE){
  if (log2.transform  != FALSE){
    dat$signal.norm <- log2(dat$signal.norm + 1)
  }
  ggplot(dat, aes(x = time, y = signal.norm, colour = experiment, group = experiment, fill = experiment)) + 
    geom_point() +
    geom_line() +
    facet_wrap(~tissue)
}



MergeRnaseqArray <- function(kallisto.long, array.long, array.var = "signal.norm"){
#   kallisto.long$experiment <- rep("rnaseq", nrow(kallisto.long))
#   array.long$experiment <- rep("array", nrow(array.long))
  
  # Get common genes
  rnaseq.genes <- unique(kallisto.long$gene)
  array.genes <- unique(array.long$gene)
  common.genes <- intersect(rnaseq.genes, array.genes)
  
  # We only want the overlap of times, RNA-Seq is subset of Array times
  rnaseq.times <- unique(kallisto.long$time)
  array.times <- unique(array.long$time)
  common.times <- intersect(rnaseq.times, array.times)
  
  kallisto.sub <- subset(kallisto.long, gene %in% common.genes & time %in% common.times)
  array.sub <- subset(array.long, gene %in% common.genes & time %in% common.times)
  
  # Fix factor levels
  kallisto.sub$gene <- factor(kallisto.sub$gene)
  array.sub$gene <- factor(array.sub$gene)
  
  # order them up
  kallisto.sub <- kallisto.sub[with(kallisto.sub, order(gene, tissue, time)), ]
  array.sub <- array.sub[with(array.sub, order(gene, tissue, time)), ]
  
  # Check my gene, tissue and times all line up
  vars <- c("gene", "tissue", "time")
  print(sapply(vars, function(var) all(kallisto.sub[[var]] == array.sub[[var]])))
  
  # Now just cbind because the gene tissue and times line up
  array.signal <- array.sub[[array.var]]
  merge.long <- cbind(kallisto.sub, array.signal)
  return(merge.long)
}

AdjustArray <- function(dat, fits.all){
  gene <- as.character(dat$gene[1])
  fit <- fits.all[gene, ]
#   if (all(is.na(fit))){
#     #     dat$signal <- NULL  # makes output DF slimmer, we only need the signal.norm column
#     #     dat$origin <- rep('unadj', nrow(dat))
#     return(data.frame(NULL))
#   }
  if (fit$model == "saturation"){
    signal.adj <- sapply(dat$signal.norm, 
                         function(m) saturation.inv2(a = fit$a, b = fit$b,k = fit$k, m))
  } else {
    signal.adj <- sapply(dat$signal.norm,
                         function(m) linear.inv2(int = fit$int, slope = fit$slope, m))
  }
  dat.new <- data.frame(gene = rep(dat$gene, 2),
                        tissue = rep(dat$tissue, 2),
                        time = rep(dat$time, 2),
                        signal.norm = c(dat$signal.norm, signal.adj),
                        origin = c(rep('unadj', nrow(dat)), rep('adj', nrow(dat))))
}

# Load matrices -----------------------------------------------------------

kallisto.path <- "data/kallisto/abundance.genecounts.matrix.txt"
array.path <- "data/array_exprs_colnames_fixed.best.probe.selected.txt"
kallisto.long <- LoadKallistoGene(kallisto.path, form = "long")
kallisto.wide <- LoadKallistoGene(kallisto.path, form = "wide")
array.long <- LoadArray(array.path, get.norm = TRUE, form = "long")
array.wide <- LoadArray(array.path, get.norm = TRUE, form = "wide")

common.samps <- intersect(colnames(kallisto.wide), colnames(array.wide))
unique.samps <- setdiff(colnames(array.wide), colnames(kallisto.wide))

# jgene <- "Elovl3"
# PlotGeneTpm(kallisto.long, jgene, log2.pseudocount=0.001, scale = 1)
# PlotGeneAcrossTissues(subset(dat.long, gene == jgene))


# # Density plot my TPM -----------------------------------------------------
# 
# scale.factor <- 100  # make kallisto range between 0 to 15 log2 scale when I do +1 pseudo count
# kallisto.long$tpm <- kallisto.long$tpm * scale.factor
# kallisto.sub <- subset(kallisto.long, tpm > 0)
# plot(density(log2(kallisto.sub$tpm + 1)))
# Get mean and variance for each gene across tissues ----------------------


# Get a model for my noise ------------------------------------------------

mean.var.gene.tiss <- array.long %>%
  group_by(gene, tissue) %>%
  summarise(Mean = mean(signal.norm),
            Var = var(signal.norm)) %>%
  .[order(.$Mean, decreasing = FALSE), ]

library(plyr)
n.per.bin <- 150
mean.var.gene.tiss$BinOrder <- factor(round_any(seq(1:nrow(mean.var.gene.tiss)), n.per.bin))
detach("package:plyr", unload=TRUE)
library(dplyr)

mean.var.meds <- mean.var.gene.tiss %>%
  group_by(BinOrder) %>%
  summarise(MedianExprs = median(Mean),
            MedianVar = median(Var))

# head(data.frame(mean.var.meds[order(mean.var.meds$MedianExprs), ]), n = 500)

fit.noise <- loess(formula = MedianVar ~ MedianExprs, 
                   data = mean.var.meds, 
                   control=loess.control(surface="direct"))

# Fit Loess Model
x <- seq(min(array.long$signal.norm), max(array.long$signal.norm), 100)  # plot full range
y <- predict(fit.noise, x)

# # Plot diagnostic plots
# pdf("plots/adjust_array_to_kallisto_diagnostics/loess_fit_diagnostics.pdf")
# plot(x, y, col='red', lwd='2', type='l', main=paste0('Loess Fit. Bin size=', n.per.bin),
#      xlab="bin.mean", ylab="bin.var",
#      xlim=c(0, 3000), ylim=c(0, 50000))
# points(mean.var.meds$MedianExprs, mean.var.meds$MedianVar)
# 
# plot(x, y, col='red', lwd='2', type='l', main=paste0('Loess Fit. Bin size=', n.per.bin),
#      xlab="bin.mean", ylab="bin.var")
# points(mean.var.meds$MedianExprs, mean.var.meds$MedianVar)
# 
# plot(x, y, col='red', lwd='2', type='l', main=paste0('Loess Fit. Bin size=', n.per.bin),
#      xlab="bin.mean", ylab="bin.var", log="xy")
# points(mean.var.meds$MedianExprs, mean.var.meds$MedianVar)
# dev.off()


# Merge RNA-Seq and array as a long data frame, separate columns ----------

merge.long <- MergeRnaseqArray(kallisto.long, array.long, array.var = "signal.norm")


fits.test <- subset(merge.long, gene %in% c(GetGeneList())) %>%
  group_by(gene) %>%
  do(FitSaturationCurve(., fit.noise, array.wide))
# 
jgene <- "Arf6"
test <- subset(merge.long, gene == "BC053393")
test <- subset(merge.long, gene == "1600029I14Rik")
test <- subset(merge.long, gene == jgene)
test.fit <- FitSaturationCurve(test, fit.noise, array.wide)

# head(merge.long)
# start <- Sys.time()  # ~4 minutes total
# fits.all <- merge.long %>%
#   group_by(gene) %>%
#   do(FitSaturationCurve(., fit.noise, array.wide))
# save(fits.all, file = "Robjs/adjust_array_to_kallisto.fits.all.Robj")
# print(Sys.time() - start)

load(file = "Robjs/adjust_array_to_kallisto.fits.all.Robj")
rownames(fits.all) <- fits.all$gene
fits.all <- data.frame(fits.all)

# # Do some diagnostic plots ------------------------------------------------
# 
gene <- "1600029I14Rik"
gene <- "4930427A07Rik"
gene <- "Arf6"

PlotDiagnostics2(gene, kallisto.wide, array.wide, fits.all)

genes <- intersect(fits.all$gene, GetGeneList())
pdf("plots/adjust_array_to_kallisto_diagnostics/diagnostics.pdf")
for (gene in genes){
  PlotDiagnostics2(gene, kallisto.wide, array.wide, fits.all)
}
dev.off()


# Adjust array to RNA-Seq -------------------------------------------------

array.adj <- array.long %>%
  subset(., gene %in% common.genes) %>%
  group_by(gene) %>%
  do(AdjustArray(., fits.all))

# Append RNA-Seq info
kallisto.long$origin <- rep("rnaseq", nrow(kallisto.long))
colnames(kallisto.long)[which(colnames(kallisto.long) == "tpm")] <- "signal.norm"  # temporarily change colnames to match array.adj

array.adj <- rbind(array.adj, kallisto.long)
colnames(kallisto.long)[which(colnames(kallisto.long) == "signal.norm")] <- "tpm"  # change it back

# Diagnostics on some genes -----------------------------------------------

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
scale_colour_manual(values=cbPalette)

array.sub <- subset(array.adj, gene %in% GetGeneList())
pdf("plots/adjust_array_to_kallisto_diagnostics/final_adj.pdf")
for (jgene in GetGeneList()){
  print(ggplot(subset(array.sub, gene == jgene), aes(x = time, y = log2(signal.norm + 0.001), colour = origin, group = origin, fill = origin)) + 
    scale_colour_manual(values=cbPalette) +
    geom_point() + geom_line() + facet_wrap(~tissue) +
    ggtitle(jgene))
}
dev.off()

pdf("plots/adjust_array_to_kallisto_diagnostics/final_adj2.pdf")
for (jgene in GetGeneList()){
  print(ggplot(subset(array.sub, gene == jgene & origin %in% c("adj", "rnaseq")), aes(x = time, y = log2(signal.norm + 0.001) , colour = origin, group = origin, fill = origin)) + 
    scale_colour_manual(values=cbPalette) +
    geom_point() + geom_line() + facet_wrap(~tissue) +
    ggtitle(jgene))
}
dev.off()


# Print out new table -----------------------------------------------------

array.adj.out <- subset(array.adj, origin == "adj") %>%
  dcast(., formula = gene ~ tissue + time, value.var = "signal.norm")

# fix colnames and set rownames
rownames(array.adj.out) <- array.adj.out$gene
array.adj.out$gene <- NULL
# Adr_24 -> Adr24 keeps life consistent
colnames(array.adj.out) <- sapply(colnames(array.adj.out), 
                                  function(s) paste(strsplit(s, "_")[[1]], collapse = ''))

# Sanity checks
rowmins <- apply(array.adj.out, 1, min)
write.table(x = array.adj.out, quote = FALSE, sep = "\t", row.names = TRUE, col.names = NA,
            file = "data/exprs_matrices/array_adj_to_kallisto.slope07.txt")
            