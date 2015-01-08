# Jake Yeung
# January 7 2015
# To eliminate genes with "noise", which shows up in projection and subsequent svd analysis, let's do projection only
# on genes that show rhythmicity in at least one tissue.
# 
# Dependencies: output of find.oscillating.genes.R is used to detect rhythmic genes.

library(ddply)
library(reshape2)
library(doMC)


# Parameters to play ------------------------------------------------------

## Specify the number of cores
mc <- 30

# Minimum log2 expression to be considered "expressed" in RNA-Seq
cutoff <- 5


clockgenes <- c('Nr1d1','Dbp', 'Arntl', 'Npas2', 'Nr1d2', 
                'Bhlhe41', 'Nfil3', 'Cdkn1a', 'Lonrf3', 
                'Tef', 'Usp2', 'Wee1', 'Dtx4', 'Asb12')

clockgenes <- c(clockgenes, 'Elovl3', 'Clock', 'Per1', 'Per2', 'Per3', 'Cry2', 'Cry1')

# My Functions ------------------------------------------------------------

scripts.dir <- "scripts"
funcs.dir <- "functions"
source(file.path(scripts.dir, funcs.dir, "MergeAndProject.R"))
source(file.path(scripts.dir, funcs.dir, "GetTissueTimes.R"))
source(file.path(scripts.dir, funcs.dir, "FourierFunctions.R"))
source(file.path(scripts.dir, funcs.dir, "LoadAndHandleData.R"))
source(file.path(scripts.dir, funcs.dir, "DataHandlingFunctions.R"))
source(file.path(scripts.dir, funcs.dir, "PlotFunctions.R"))
source(file.path(scripts.dir, funcs.dir, "GetTopNValues.R"))
source(file.path(scripts.dir, funcs.dir, "OuterComplex.R"))
source(file.path(scripts.dir, funcs.dir, "OrderPhaseMatrix.R"))

ProjectToFrequency <- function(df, my.omega, normalize = TRUE){
  # Perform fourier transform and normalize across all frequencies.
  # 
  # Input:
  # df: long dataframe containing expression and time columns
  # to be transformed to frequency domain.
  # 
  # omega in which we are interested.
  
  # if we normalize, get list of omegas.
  # otherwise we just use omega
  if (normalize){
    # t <- sort(unique(df$time))
    # n.timepoints <- length(t)
    # interval <- t[2] - t[1]
    omegas <- GetOmegas(remove.zero = TRUE)
  } else {
    omegas <- my.omega
  }
  
  sum.transforms <- 0
  
  for (omega in omegas){
    # take gene from a tissue and get its fourier transform
    transformed <- df$exprs %*% exp(-1i * omega * df$time)
    # normalize by number of datapoints
    transformed <- transformed / length(df$exprs)
    sum.transforms <- sum.transforms + Mod(transformed)
    
    if (omega == my.omega){
      my.transformed <- transformed
    }
  }
  
  if (normalize){
    # Normalize across omegas
    # Multiply by mean
#     jmean <- mean(subset(df, experiment == "rnaseq")$exprs)  # rnaseq exprs more reliable than microarray.
#     factor <- jmean
    
    # if median is 0, then set factor to 0, otherwise 1
    factor <- 1
    cutoff <- 5
    jmedian <- median(subset(df, experiment == "rnaseq")$exprs)
    if (jmedian <= cutoff){
      factor <- 0
    }
    
    my.transformed <- (my.transformed / sum.transforms) * factor
  } 
  
  data.frame(exprs.transformed = my.transformed)
}

DoFourier <- function(exprs, time, omega){
  p <- exprs %*% exp(-1i * omega * time)
  return(p)
}

GetOmegas <- function(n.timepoints=24, interval=2, remove.zero=FALSE){
  # Get omegas from a time series
  # 
  # remove.zero: removes first element (omega = 0).
  # Useful for normalizing Fourier when I don't want omega = 0 case
  
  # Output:
  # list of omegas: can be used for fourier transform
  
  t <- seq(from = 0, by = 1, length.out = n.timepoints)  # will account for interval later
  
  # get harmonic frequencies
  # this creates harmonic frequencies from 0 to .5 in steps of 1/n.timepoints
  
  mid.point <- length(t) / 2 + 1
  freqs <- t[1:mid.point] / n.timepoints
  
  # Convert to period before adjusting for interval
  T.unadj <- 1 / freqs
  
  # Account for interval
  T <- T.unadj * interval
  
  # Calculate omegas
  omegas <- 2 * pi / T
  
  if (remove.zero){
    omegas <- omegas[2:length(omegas)]
  }
  
  return(omegas)
}

ConvertLongToWide <- function(long.df, measurement.var = "exprs.transformed"){
  wide.df <- dcast(dat.proj, gene ~ tissue, value.var = measurement.var)
  # first row is gene name, let's make them rowname and remove first column.
  
  rownames(wide.df) <- wide.df$gene
  
  wide.df <- subset(wide.df, select = -c(gene))
  
  return(wide.df)
}

PlotGeneAcrossTissues <- function(dat, jtitle){
  m <- ggplot(dat, aes(x = time, y = exprs,
                       group = experiment, 
                       colour = experiment)) + 
    geom_point() + 
    geom_line() + 
    facet_wrap(~tissue) +
    ggtitle(jtitle)
  return(m)
}



# Define dirs -------------------------------------------------------------

# define dirs
data.dir <- "data"
normalized.array.fname <- "array.adj.0.07.txt"
normalized.array.path <- file.path(data.dir, normalized.array.fname)
rna.seq.fname <- "rna_seq_deseq_counts_colnames_fixed.txt"
rna.seq.path <- file.path(data.dir, rna.seq.fname)

# Load file ---------------------------------------------------------------

normalized.array <- LoadNormalizedArray(normalized.array.path)
rna.seq.exprs <- LoadRnaSeq(rna.seq.path)


# Log2 transform of array and rnaseq --------------------------------------

normalized.array <- log2(normalized.array + 1)
rna.seq.exprs <- log2(rna.seq.exprs + 1)


# Remove genes from RNA-Seqs that are not expressing. ---------------------

# plot and show cutoff visually
plot(density(unlist(rna.seq.exprs)))
abline(v = cutoff, col = 'red')

jmean <- apply(rna.seq.exprs, 1, mean)

genes.not.expressed <- names(jmean[which(jmean < cutoff)])
all.genes <- rownames(rna.seq.exprs)
genes.expressed <- all.genes[which(!all.genes %in% genes.not.expressed)]

print(paste0(length(genes.not.expressed), 
             " not expressed at cutoff = mean(", cutoff, 
             ") for a gene across all tissues. Removing..."))

# Remove unexpressed genes

rna.seq.exprs.filtered <- rna.seq.exprs[genes.expressed, ]


# Take only common genes --------------------------------------------------

array.genes <- rownames(normalized.array)
rna.seq.genes <- rownames(rna.seq.exprs.filtered)

common.genes <- intersect(array.genes, rna.seq.genes)

print(paste(length(common.genes), "common genes between array and rnaseq"))

normalized.array <- normalized.array[common.genes, ]
rna.seq.exprs.filtered <- rna.seq.exprs.filtered[common.genes, ]

# Merge data into long format ---------------------------------------------

dat <- MergeToLong(normalized.array, rna.seq.exprs.filtered)

dat.split <- split(dat, dat$tissue)


# Project my data ---------------------------------------------------------

omega <- 2 * pi / 24
omegas <- GetOmegas()

# # # test my function
# clockgenes <- c(clockgenes, "Spint4", "Defb23", "Defb34", "Gje1", "Gm10233", "Mir181b.2", "n.R5s54")
# problem genes
problem.genes <- c("Hephl1", "Dnmt3l", "Gm10804", "Fgf14", "Spint4", "Defb23", "Defb34", "B3galt2")
problem.genes <- c("B3galt2", "Sgcg", "Rgs16", "Ddit4", "Smpx", "Eef1a2", "D3Bwg0562e", "Slc9a2", 
                   "Tceal3", "Murc", "Fgf10", "Dclk1", "Ccne1", "Cox6a2", "Trim63", 
                   "Kcng4", "Slc17a9", "Diras2", "Txlnb", "A830018L16Rik", "Elovl3", "Alb")
problem.genes <- c("Diras2", "Kcng4", "Trim63", "Cox6a2")
df <- subset(dat.split$Liver, gene %in% c(clockgenes, problem.genes))
df.proj <- ddply(df, .(gene), ProjectToFrequency, omega, normalize = TRUE)
df.proj$mod <- Mod(df.proj$exprs.transformed)
df.proj <- df.proj[order(df.proj$mod, decreasing = TRUE), ]
(df.proj)
# 
# # test individual genes
# jgene <- "Gm10233"
# jgene <- "Mir181b.2"
# jgene <- "Gje1"
jgene <- "B3galt2"
jgene <- "Arntl"
jmean <- mean(subset(dat.split$Liver, gene == jgene & experiment == "rnaseq")$exprs)
sapply(omegas, DoFourier, exprs = subset(dat.split$Liver, gene == jgene)$exprs, time = subset(dat.split$Liver, gene == jgene)$time)


start.time <- Sys.time()
if (getDoParWorkers() == 1){
  registerDoMC(mc)
}
print(paste("Parallel processing on:", getDoParWorkers(), "cores"))
dat.split.proj <- lapply(dat.split, function(x){
  ddply(x, .(gene), ProjectToFrequency, my.omega = omega, normalize = TRUE, .parallel = TRUE)
})
print(Sys.time() - start.time)


# Add tissue information into each list -----------------------------------

for (tissue in names(dat.split.proj)){
  dat.split.proj[[tissue]]$tissue <- tissue
}

# Combine data ------------------------------------------------------------

dat.proj <- do.call(rbind, dat.split.proj)


# Remove NaNs -------------------------------------------------------------

# dat.proj <- dat.proj[is.finite(dat.proj$exprs.transformed), ]

# Clean up ----------------------------------------------------------------

# rm(dat.split.proj, dat.split)

# Do SVD on this normalized dataset ---------------------------------------

# long to wide conversion
dat.wide <- ConvertLongToWide(dat.proj, measurement.var = "exprs.transformed")

# Complete cases only
dat.wide <- dat.wide[complete.cases(dat.wide), ]

s <- svd(dat.wide)


# Set rownames for u and v ------------------------------------------------

rownames(s$u) <- rownames(dat.wide)
rownames(s$v) <- colnames(dat.wide)

# Plot interesting components ---------------------------------------------

sing.vals <- seq(length(s$d))
# sing.vals <- c(1, 2, 3)

# ScreePlot
plot(s$d^2 / sum(s$d ^ 2), type='o')  # Manual screeplot. 

for (sing.val in sing.vals){
  print(paste("Plotting component:", sing.val))
  eigengene <- s$v[, sing.val]
  jmax <- max(Mod(eigengene))
  PlotComplex(eigengene, 
              axis.min = -jmax,
              axis.max = jmax,
              labels = tissues, 
              col = "HSV",
              add.text.plot = FALSE, 
              main = paste("Component:", sing.val),
              jpch = 20,
              threshold = 0)
  
  eigensample <- s$u[, sing.val]
  jmax <- max(Mod(eigensample))
  PlotComplex(eigensample,
              axis.min = -jmax,
              axis.max = jmax,
              labels = names(eigensample),
              col = "HSV",
              main = paste("Component:", sing.val),
              add.text.plot = FALSE,
              jpch = 1,
              threshold = 0.5 * jmax,
              verbose = FALSE)
  
  # show how fast the values drop
  top.genes.all <- GetTopNValues(Mod(eigensample), N = length(eigensample))
  plot(top.genes.all$vals, 
       log = "x", 
       main = paste("Component:", sing.val),
       ylab = "Modulus of gene loading",
       xlab = "Gene index")
  # Plot args in color
  top.genes <- GetTopNValues(Mod(eigensample), N = 100)# list of $vals $i
  # use drop to keep rownames
  outer.prod.mat <- s$d[sing.val] * OuterComplex(s$u[top.genes$i, sing.val, drop = FALSE], t(s$v[, sing.val, drop = FALSE]))
  outer.prod.mat <- OrderPhaseMatrix(outer.prod.mat)
  PlotArgsMatrix(outer.prod.mat, main = paste("Component", sing.val))
  
  # print top 100 genes (copy and paste-able)
  cat(paste0(head(names(top.genes$vals), n = 50), "\n"))
  
#   # print top 20 genes (list-able)
  cat(paste(head(names(top.genes$vals), n = 20), collapse = '", "'))
  cat("\n")
}


# Plot interesting genes --------------------------------------------------

# component 2: normalized and multiplied and filtered (these are LIVER GENES)
genes <- c("Elovl3", "Ddit4", "Chka", "Pdk4", "Slc45a3", "Bhmt", "Cml5", "Osgin1", "Ppard", "BC029214", "Rgs16", "Acacb", "Leap2", "Upp2", "Ethe1")

genes <- c("Cmah", "Cyp2f2", "Nr1d1", "C4b", "Slc25a10", "Cd37", "Dbp", "Zbtb16", "Ttr", "Rbp4", "Aplnr", 
           "Pnpla3", "Serpina3n", "Apoa2", "Sucnr1", "Alb", "Cfb", "Sncg", "Aldh1a1", "Cyp2b10")

# component 2: normalized, no multiply, filter by MEAN 
genes <- c("B3galt2", "Sgcg", "Rgs16", "Ddit4", "Smpx", "Eef1a2", "D3Bwg0562e", "Slc9a2", "Tceal3", "Murc", "Fgf10", "Dclk1")

# component 2: normalized. multiply by 1 or 0, filter by mean
genes <- c("Rgs16", "Ddit4", "Tceal3", "Dclk1", "Ccne1", "Diras2", "Slc17a9", 
           "Cntfr", "Lrrc2", "Txlnb", "Aqp7", "Ppard", "Mb", "Lrfn3", "Acot1", 
           "Kif26b", "Dtna", "Slc41a3", "Ebf1", "Pdk4")

# component 2: normalized, multiply by 1 or 0 (if less than cutoff) filter by mean
genes <- c("Rgs16", "Ddit4", "Ppard", "Cntfr", "Lonrf3", "Lrfn3", "Slc17a9", "Slc45a3", "Chka", "Pdk4", "Phospho1", 
           "BC029214", "Osgin1", "I830012O16Rik", "Mpzl1", "Celsr1", "Arsg", "Abtb2", "Lgalsl", "Fbxo44")

genes <- c("Npas2", "Cmah", "Mthfd1l", "Dbp", "Cfb", "Agt", "Lonrf3", "Nr1d1", "Usp13", "Eno3", "Cldn1", 
           "Slc38a3", "Serpina3n", "Hif3a", "Rasl11a", "Sncg", "Serpina3c", "Slc25a10", "Cyp2f2", "Tcea3")
dat.sub <- subset(dat, gene %in% genes)
dat.sub$gene <- factor(dat.sub$gene)

for (jgene in genes){
  print(jgene)
  print(PlotGeneAcrossTissues(subset(dat.sub, gene == jgene), jtitle = jgene))
}

# # Interesting genes: Fam150b, Cml5, Tshr
# 
# # genes <- c("Dbp", "Defb22", "Defb23", "Defb34", "Npas2", "Nr1d1", "Olfr1014", "Olfr1015", "Olfr1052", "Olfr731", "Spint4")
# 
# genes <- c("Gje1", "Gm10233", "Lonrf3", "Mir181b.2", "n.R5s54")
# 
# genes <- c("Cry2", "Arntl", "Dbp", "Gje1", "n.R5s54", "Npas2")
# 
# genes <- c("Abcd2", "Abcg5", "Abcg8", "Abtb2", "Acacb", "Acnat2", "Acot1", "Acsl1", "Adck3", "Adprhl1", 
#            "Alas1", "Amer1", "Angptl4", "Aqp3", "Aqp8", "Arrdc3", "Arsg", "Avpr1a", "B3galt2", "BC029214", 
#            "Bhmt", "Camk2a", "Ccrn4l", "Celsr1", "Chka", "Ckm", "Clec2h", "Clpx", "Cml5", "Cntfr", 
#            "Col11a1", "Col27a1", "Cox6a2", "Crip2", "Ddc", "Ddit4", "Dupd1", "Ebf1", "Eef1a2", "Elmod3", 
#            "Elovl3", "Ethe1", "Fdft1", "Fgf14", "Gm10762", "Gm10804", "Gm11437", "Gm15998", "Gm4952", 
#            "Gpam", "Gpcpd1", "Gprin3", "Hacl1", "Hdhd3", "Hsd3b7", "I830012O16Rik", "Inca1", "Insc", "Insig2", "Irf6", "Jph3", "Klb", "Klhl30", "Leap2", "Lep", "Lgalsl", "Loxl4", "Lrfn3", "Lrrc30", "Mb", "Myc", "Myh6", "Myh8", "Ndrg1", "Net1", "Nr1d1", "Nrk", "Osgin1", "Paqr9", "Pck1", "Pdk4", "Pfkfb3", "Plin5", "Plxna2", "Pnp", "Pnp2", "Ppard", "Pxmp4", "Rassf6", "Rgs16", "Sco2", "Sept9", "Slc17a9", "Slc2a2", "Slc38a2", "Slc38a3", "Slc45a3", "Slc7a2", "Smagp", "Smpx", "Smtnl2", "Snhg3", "Spata22", "St5", "St6galnac6", "Tbx15", "Tfpi2", "Tfrc", "Tjp3", "Tmem37", "Tmie", "Tnfrsf19", "Tra2a", "Trdn", "Trim63", "Tshr", "Tspan4", "Tst", "Tubb2a", "Ube2u", "Upp2", "X1300002K09Rik", "X4930444P10Rik", "Xirp2")
# 
# genes <- c("Ddit4", "Elovl3", "Fgf14", "Pdk4", "Tshr")  # normalizing and multiplying by mean
# genes <- c("Fam124b", "Fam150b", "Lrrc30", "Myh8", "Tshr", "X4930444P10Rik")  # removing "unexpressed genes" the normalizing, no multipy by mean
# genes <- c("Tshr", "Fam124b", "Myh8", "Lrrc30", "X4930444P10Rik", "Fam150b", "B3galt2", "Dnmt3l", 
#            "Hephl1", "Gm10804", "Tacr3", "Krt9", "Tmem72", "CN725425", "Tbx15", "Adprhl1", 
#            "Cyp11a1", "Cpa4", "Ripply1", "Tpo", "Mc5r", "Ano7", "Wnt16", "Rab9b", "Lsmem1", 
#            "Bcl2l10", "Nr1d1", "Cml5", "Ccdc121", "Kcnq3", "Sgcg", "Lonrf3", "Aqp3", "Cd209d", 
#            "Myf5", "Pcdh8", "Kcne1l", "Cyp24a1", "Casr", "Tpsb2", "Smpx", "Dmrtc1a", "Defb2")
# 
# # 2nd component genes: normalize and multiply by mean
# genes <- c("Bhmt", "Chka", "Ddit4", "Elovl3", "Fgf14", "Pdk4", "Rgs16", "Tshr") 
# 
# # 2nd component genes: normalize and multiply by mean (RNASEQ ONLY)
# genes <- c("Elovl3", "Ddit4", "Chka", "Pdk4", "Slc45a3", "Bhmt", "Cml5", 
#            "Osgin1", "Ppard", "BC029214", "Rgs16", "Acacb", 
#            "Leap2", "Upp2", "Ethe1", "Slc30a10", "Tst", "Aqp8", 
#            "Tshr", "Celsr1")
# 
# # 3rd component genes: normalize and multiply by mean of RNASEQ
# genes <- c("Adrb3", "Ahsg", "Alb", "Aldh1a1", "Aplnr", "Apoa2", 
#            "Bhlhe40", "Bhlhe41", "C4b", "Cd37", "Cfb", "Cmah", 
#            "Cml5", "Cxcr4", "Cyp2b10", "Cyp2f2", "Dbp", "Dgat2", 
#            "Ebf1", "Eno3", "Fbxo21", "Gsta3", "Hdc", "Hspa5", 
#            "Igfbp4", "Lrg1", "Nnat", "Npas2", "Nr1d1", "Pnpla3", 
#            "Rbp4", "Retn", "Serpina3n", "Slc25a1", "Slc25a10", 
#            "Sncg", "Sucnr1", "Trf", "Ttr", "Zbtb16")
# 
# # problem genes
# genes <- c(c("Hephl1", "Dnmt3l", "Gm10804", "Fgf14"))
# genes <- c("Spint4", "Defb23", "Defb34")
# 
# # 2nd component: no normalization (NOISY GENES)
# genes <- c("Adam7", "Arntl", "Cst8", "Dbp", "Defb47", "Defb48", "Gpx5", "Lcn10", "Lcn8", "Lcn9", "Ly6g5b", "Ly6g5c", "Npas2", "Spint4", "Svs2")
# 
# # 2nd component: normalization only, no mean multiplication
# genes <- c("Fam124b", "Fam150b", "Lrrc30", "Myh8", "Tshr", "X4930444P10Rik")
# 
# # 2nd component: no division, just multiply by mean
# genes <- c("Adam7", "Bhmt", "Cst11", "Cst8", "Defb20", "Defb48", "Elovl3", 
#            "Gpx5", "Lcn8", "Lcn9", "Ly6g5b", "Rnase10", "Svs2", "X5830403L16Rik", "X9230104L09Rik")
# 
# dat.sub <- subset(dat, gene %in% genes)
# for (jgene in genes){
#   print(jgene)
#   m <- PlotGeneAcrossTissues(subset(dat.sub, gene == jgene), jtitle = jgene)
#   print(m)
# }


