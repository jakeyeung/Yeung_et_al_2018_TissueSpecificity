# Jake Yeung
# find oscillating genes
# find.oscillating.genes.R
# Dec 4 2014

library(ggplot2)
library(plyr)
library(foreach)  # for parallelization
library(doParallel)  # for parallelization

# Adjustable parameters ---------------------------------------------------

outfile <- "data/oscillating_genes.txt"
gene.header <- "gene"
params.header <- c("amp", "phase", "int.array", "int.rnaseq")

# Define constants --------------------------------------------------------

pval.threshold <- 1e-6
amp.threshold <- 1.5
T <- 24  # hours
omega <- 2 * pi / T


# Define genes ------------------------------------------------------------

clockgenes <- c("Nr1d1","Dbp","Arntl","Npas2","Nr1d2","Bhlhe41","Nfil3",
                "Cdkn1a","Lonrf3","Tef","Usp2","Wee1","Dtx4","Asb12",
                "Elovl3","Clock","Per1","Per2","Per3","Cry2","Cry1")
tissuegenes <- c("Plbd1","Acacb","Hadh","Decr1","Acadl","Ech1","Acsl1","Cpt2","Acaa2")  # PCA 1 from explore_data.R
tissuegenes <- c(tissuegenes, "Elovl1","Slc35b3","Fkbp15","Slc39a8","Sep15","Mospd2","Med11","Slco2a1","Arf6","Yipf6")  # PCA 2 from explore_data.R



# Functions ---------------------------------------------------------------

scripts.dir <- 'scripts'
funcs.dir <- 'functions'
source(file.path(scripts.dir, funcs.dir, "ReplaceNegs.R"))
source(file.path(scripts.dir, funcs.dir, "GetTissueTimes.R"))  # get tissue and times
source(file.path(scripts.dir, funcs.dir, "DataHandlingFunctions.R"))
source(file.path(scripts.dir, funcs.dir, "RegressionFunctions.R"))

# Fit complex model
FitComplexModel <- function(df, omega = 2 * pi / 24){
  lm(exprs ~ 0 + experiment + cos(omega * time) + sin(omega * time), 
     data = df)
}

# Get parameters from complex model
GetParamsComplexModel <- function(myfit){
  model.params <- coef(myfit)
  return(list(intercept.array = model.params[1],
              intercept.rnaseq = model.params[2],
              a = model.params[3],
              b = model.params[4]))
}

# Fit flat model
FitFlatModel <- function(df){
  fit.flat <- lm(exprs ~ 0 + experiment, data =df)
}

# Get pval
GetFtestPval <- function(myfit){
  pval <- myfit["Pr(>F)"][[1]][2]
  return(pval)
}

# find rhythmic
FindRhythmic <- function(df){
  fit.complex <- FitComplexModel(df)
  fit.flat <- FitFlatModel(df)
  ftest.flat.or.rhythmic <- anova(fit.flat, fit.complex)
  pval <- GetFtestPval(ftest.flat.or.rhythmic)
  return(list(pval=pval, fit=fit.complex))
}

# get amplitude and phase
GetAmpPhase <- function(df){
  fit.complex <- FitComplexModel(df)
  fit.flat <- FitFlatModel(df)
  ftest.flat.or.rhythmic <- anova(fit.flat, fit.complex)
  pval <- GetFtestPval(ftest.flat.or.rhythmic)
  params <- GetParamsComplexModel(fit.complex)  # model = intercept.array + intercept.rnaseq + acos(t) + bsin(t)
  # for amp and phase calculation: see http://www.intmath.com/analytic-trigonometry/6-express-sin-sum-angles.php
  amp <- unname(sqrt(params[["a"]] ^ 2 + params[["b"]] ^ 2))
  phase <- unname(atan2(params[["a"]], params[["b"]]))
  # phase <- atan(params[["b"]] / params[["a"]])
  int.array <- params[["intercept.array"]]
  int.rnaseq <- params[["intercept.rnaseq"]]
  return(list(pval=pval,
         amp=amp, 
         phase=phase, 
         int.array=int.array, 
         int.rnaseq=int.rnaseq))
}

# list to textfile
ListToTextfile <- function(mylist, myfile, ncols=5){
  lapply(mylist, write, myfile, append=TRUE, ncolumns=ncols)
}

# oscillating function
oscillate <- function(params, omega, t){
  params[1] + params[2] * cos(omega * t) + params[3] * sin(omega * t)
}

# guess if simple or complex from fit parameters
IsSimple <- function(fit, simple.df=3, complex.df=4){
  # If fit has simple.df parameters, it is simple. Return TRUE
  # if fit has complex.df parameter,s it is complex. Return FALSE
  # 
  # Args:
  # fit from lm
  # 
  # Return:
  # TRUE if simple
  # FALSE if complex
  # NA if neither, print warning.
  n.params <- summary(fit)$df[[1]]
  if (n.params == simple.df){
    return(TRUE)
  } else if (n.params == complex.df){
    return(FALSE)
  }
    else 
      warning(paste("Fit parameters neither has", simple.df, 
                    "or", complex.df, "parameters. 
                    It has", n.params, ". Defaulting to FALSE."))
  return(FALSE)
}

GetHeader <- function(gene.header, params.header, tissues){
  # Create column names for each parameter in complex fit for each tissue
  # First element in header is gene.header.
  
  # Make header like Tissue:Parameter
  params.tissue.header <- paste0(rep(tissues, each = length(params.header)), ":", params.header)
  
  # Add gene to it and return
  header <- c(gene.header, params.tissue.header)
  return(header)
}

IsRhythmic <- function(fit.params, pval.threshold, amp.threshold=0){
  # Given parameters from complex fit (we know pval and amp) return TRUE
  # if it passes threshold. Return FALSE otherwise.
  pval <- fit.params[["pval"]]
  amp <- fit.params[["amp"]]
  if (pval <= pval.threshold & amp >= amp.threshold){
    return(TRUE)
  } 
  return (FALSE)
}

GetRhythmicGene <- function(fit.params, pval.threshold, amp.threshold){
  if (IsRhythmic(fit.params, pval.threshold, amp.threshold)){
    return(TRUE)
  }
  return(FALSE)
}

GetFitInfo <- function(fit.list.tissue.gene, param){
  return(fit.list.tissue.gene[[param]])
}

GetParamsVector <- function(fit.list.tissue, param="amp"){
  # Get parameters such as amp, phase, pval, int.array, int.rnaseq
  # across all tissues.
  # 
  # param = "amp" | "pval" | "phase" | "int.array" | "int.rnaseq"
  #
  # Returns list of lists, one of each tissue
  params.vector <- lapply(fit.list.tissue, GetFitInfo, param)
  return(unlist(params.vector))
}

PlotAcrossTissues <- function(dat, fit.list, jgene, jtitle){
  # Plot expression across tissues for a gene.
  # TODO: Find a way to incorporate fit.list parameters into the data
  
  if (missing(jtitle){
    jtitle <- jgene
  })
  
  dat.sub <- subset(dat, gene == jgene)
  
  # Create new labels
  labels <- vector(mode = "list", length = length(levels(dat.sub$tissue)))  # init
  names(labels) <- levels(dat.sub$tissue)
  
  for (tissue in levels(dat.sub$tissue)){
    params <- fit.list[[tissue]][[jgene]]
    pval <- signif(params[["pval"]], 2)
    amp <- signif(params[["amp"]], 2)
    phase <- signif(params[["phase"]], 2)
    int.array <- signif(params[["int.array"]], 2)
    int.rnaseq <- signif(params[["int.rnaseq"]], 2)
    
    label <- paste0(c(tissue, 
                      paste0("pval=", pval), 
                      "\n",
                      paste0("amp=", amp),
                      "\n",
                      paste0("phase=", phase)),
                      collapse = ",")
    
    labels[[tissue]] <- label
  }
  
  tissue.labels <- character(length(dat.sub$tissue))
  i <- 1
  for (tissue.label in dat.sub$tissue){
    tissue.labels[i] <- labels[[dat.sub$tissue[i]]]
    i <- i + 1
  }
  
  dat.sub$label <- tissue.labels
  
  m <- ggplot(dat.sub, aes(x = time, y = exprs, 
                           group = experiment, 
                           colour = experiment))
  m <- m + geom_point() + 
    geom_line() + 
    facet_wrap(~label) +
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


# Define plot files -------------------------------------------------------

fit.plot <- "plots/fourier.analysis.top.oscillators.pdf"

# Load file ---------------------------------------------------------------

normalized.array <- read.table(normalized.array.path)
rna.seq.exprs <- read.table(rna.seq.path, header=TRUE, sep='\t')

# Handle duplicate rownames: RNASEQ ------------------------------------

rownames(rna.seq.exprs) <- make.names(rna.seq.exprs$gene, unique=TRUE)

drop.cols <- c("gene")
rna.seq.exprs <- rna.seq.exprs[, !(names(rna.seq.exprs) %in% drop.cols)]

Peek(rna.seq.exprs)  # expect gene names as row names, tissues in columns



# How many have negative values? ------------------------------------------

negs <- apply(normalized.array, 1, function(x){
  if (min(x) < 0){
    return(1)
  } else {
    return(0)
  }
})

problem.genes <- names(negs[which(negs == 1)])

rm(negs)


# Remove problem genes from analysis --------------------------------------

genes <- rownames(normalized.array)
filtered.genes <- genes[which(!genes %in% problem.genes)]
normalized.array <- normalized.array[filtered.genes, ]
rna.seq.exprs <- rna.seq.exprs[filtered.genes, ]

# Log2 transform of array and rnaseq --------------------------------------

normalized.array <- log2(normalized.array + 1)
rna.seq.exprs <- log2(rna.seq.exprs + 1)


# Merge data into long format ---------------------------------------------

tissues <- GetTissues(colnames(normalized.array))
times.array <- GetTimes(colnames(normalized.array))
times.rnaseq <- GetTimes(colnames(rna.seq.exprs))


long.array <- data.frame(gene=rep(filtered.genes, length(tissues) * length(times.array)),
                         tissue=as.factor(rep(tissues, each=(length(filtered.genes) * length(times.array)))),
                         time=as.numeric(rep(times.array, each=length(filtered.genes))),
                         experiment=as.factor("array"),
                         exprs=unlist(normalized.array))

long.rnaseq <- data.frame(gene=rep(filtered.genes, length(tissues) * length(times.rnaseq)),
                          tissue=rep(tissues, each=(length(filtered.genes) * length(times.rnaseq))),
                          time=as.numeric(rep(times.rnaseq, each=length(filtered.genes))),
                          experiment=as.factor("rnaseq"),
                          exprs=unlist(rna.seq.exprs))

dat <- rbind(long.array, long.rnaseq)

str(dat)

# cleanup
rm(long.array, long.rnaseq)


# Do my lm fit ------------------------------------------------------------
# slow: about 8 minutes

#setup parallel backend to use 12 processors
cl<-makeCluster(12)
registerDoParallel(cl)

print('Begin fitting lm models')
print(Sys.time())
fit.list <- foreach(i = 1:length(tissues), 
                     .packages="plyr") %dopar% {
                       jtissue <- tissues[i]
                       dat.tiss <- subset(dat, tissue %in% c(jtissue))
                       fit.out <- dlply(dat.tiss, .(gene), GetAmpPhase)
                     }
names(fit.list) <- tissues
print(Sys.time())
print('End fitting lm models')

# Write list to textfile -- ------------------------------------------------

print(paste("Writing fit results to file.", Sys.time()))
# Write header
myheader <- GetHeader(gene.header, params.header, tissues)
write(myheader, outfile, ncolumns = length(myheader), append = TRUE, sep = "\t")

# slow, takes about 2 minutes
for (gene in filtered.genes){
  # init vector
  writerow <- character(length(myheader))
  writerow[1] <- gene  # first column is gene
  index.start <- 2  # begin tissue parameters at column 2
  index.jump <- length(params.header) + 1  # step after each loop to write to next part of vector
  index.end <- index.start + length(params.header)
  for (tissue in tissues){
    writerow.temp <- unlist(fit.list[[tissue]][[gene]], use.names = FALSE)
    # Append to vector
    writerow[index.start:index.end] <- writerow.temp
    index.start <- index.start + index.jump
    index.end <- index.end + index.jump
  }
  # write to file
  write(writerow, outfile, ncolumns = length(myheader), append = TRUE, sep = "\t")
}
rm(writerow.temp)
print(paste("Done writing fit results to file.", Sys.time()))

# Create data frame for ggplot --------------------------------------------

fits.df <- data.frame(tissue=character(),
                      gene=character(),
                      amp=character(),
                      phase=character(),
                      pval=character())

params <- c("amp", "phase", "pval")
for (tissue in tissues){
  df.temp <- data.frame(tissue = rep(tissue, length(filtered.genes)), 
                        gene = names(fit.list[[tissue]]),
                        amp = vector(length = length(filtered.genes)),
                        phase = vector(length = length(filtered.genes)),
                        pval = vector(length = length(filtered.genes)))
  
  for (jparam in params){
    params.vector <- GetParamsVector(fit.list[[tissue]], param = jparam)
    # append to df
    df.temp[[jparam]] <- params.vector
  }
  # row append
  fits.df <- rbind(fits.df, df.temp)
}
rm(df.temp)


# Plot amplitude distributions --------------------------------------------

pdf('plots/amplitude_distributions_across_tissues.pdf')
ggplot(subset(fits.df, pval <= pval.threshold), aes(x = amp, fill = tissue)) + 
  geom_density() + 
  facet_wrap(~tissue)
dev.off()


# Plot phase distributions ------------------------------------------------

pdf('plots/phase_distributions_across_tissues.pdf')
ggplot(subset(fits.df, pval <= pval.threshold), aes(x = phase, fill = tissue)) + 
  geom_density() + 
  facet_wrap(~tissue)
dev.off()



# Get top rhythmic genes --------------------------------------------------

fits.df.subset <- subset(fits.df, pval < pval.threshold)

# order by tissue, then by decreasing amp
fits.df.subset <- fits.df.subset[order(fits.df.subset$tissue, 
                                       -fits.df.subset$amp), ]

# get median amp for each tissue
amp.summary <- ddply(fits.df.subset, .(tissue), summarise, 
                     pval.min = min(pval),
                     pval.max = max(pval),
                     pval.med = median(pval),
                     amp.min = min(amp),
                     amp.max = max(amp),
                     amp.med = median(amp),
                     amp.quantile = quantile(amp, prob = 0.75))


fits.df.subset.top.n <- data.frame(tissue = character(),
                                   gene = character(),
                                   amp = character(),
                                   phase = character(),
                                   pval = character())

for (jtissue in tissues){
  amp.quantile <- amp.summary[amp.summary$tissue == jtissue, "amp.quantile"]
  tissue.temp <- subset(fits.df.subset, tissue == jtissue & amp >= amp.quantile)
  fits.df.subset.top.n <- rbind(fits.df.subset.top.n, tissue.temp)
}

head(fits.df.subset.top.n)
tail(fits.df.subset.top.n)


# Plot amplitudes ---------------------------------------------------------

ggplot(subset(fits.df.subset.top.n, pval <= pval.threshold), aes(x = amp, fill = tissue)) + 
  geom_density() + 
  facet_wrap(~tissue)


# Summarize by tissue and genes -------------------------------------------


# summmarise by tissue
tissues.sum <- ddply(fits.df.subset.top.n, .(tissue), summarise, 
                    count = length(amp), 
                    avg.amp = mean(amp), 
                    min.amp = min(amp),
                    max.amp = max(amp),
                    min.pval = min(pval),
                    max.pval = max(pval))
(tissues.sum)

# summarise by gene
r.genes.sum <- ddply(fits.df.subset.top.n, 
                     .(gene), 
                     summarise,
                     count = length(amp),
                     avg.amp = mean(amp),
                     min.amp = min(amp),
                     max.amp = max(amp),
                     min.pval = min(pval),
                     max.pval = max(pval))

# order by counts (descending), followed by min.pval (ascending), followed by max amp (descending)
r.genes.sum <- r.genes.sum[order(-r.genes.sum$count, r.genes.sum$min.pval, -r.genes.sum$max.amp), ]

str(r.genes.sum)
head(r.genes.sum)


# Plot distribution of counts ---------------------------------------------

pdf("plots/distribution_of_tissue_specificness_rhythmic_genes.pdf")
ggplot(r.genes.sum, aes(x = count)) + 
  geom_histogram() + 
  scale_x_discrete(name="Number of tissues showing rhythmicity for a gene") + 
  scale_y_continuous(name="Number of genes")
dev.off()


# Plot distribution of tissues for counts == 1 ----------------------------

# peripheral genes: genes where rhythmic only in one tissue but not others.
peripheral.genes <- subset(r.genes.sum, count == 1)$gene

subset.periph <- subset(fits.df.subset.top.n, gene %in% peripheral.genes)

subset.periph.summ <- ddply(subset.periph, .(tissue), summarise, count = length(gene))

subset.periph.summ <- subset.periph.summ[order(subset.periph.summ$count, decreasing = TRUE), ]

(subset.periph.summ$tissue)

# relevel tissues by decreasing counts

subset.periph$tissue <- factor(subset.periph$tissue, levels = subset.periph.summ$tissue)

# now plot

pdf("plots/number_of_tissue_specific_rhythmic_genes.pdf")
ggplot(subset.periph, aes(x = tissue)) + geom_bar(stat = "bin") + ggtitle("Number of tissue-specific rhythmic genes across tissues")
dev.off()

# Plot selected genes -----------------------------------------------------

my.genes <- r.genes.sum$gene

dat.sub.temp <- subset(dat, gene %in% my.genes)  # for speed?

print(paste('Printing top genes.', Sys.time()))
pdf("plots/rhythmic_genes_across_tissues.pdf")
# slow takes about 15 minutes depending on how many genes you got...
gene.count <- 0
for (gene in my.genes){
  if (gene.count %% 100 == 0){
    print(paste0(gene.count, "/", length(my.genes), " ", Sys.time()))
  }
  gene.count <- gene.count + 1
  
  count <- r.genes.sum[r.genes.sum$gene == gene, "count"]
  jtitle <- paste(gene, "rhythmic in", count, "tissues")
  
  m <- PlotAcrossTissues(dat.sub.temp, fit.list, gene, jtitle)
  print(m)
}
dev.off()
print(paste('Done printing top genes.', Sys.time()))

rm(dat.sub.temp)

# parallelization fail

# #setup parallel backend to use 30 processors
# cl<-makeCluster(30)
# registerDoParallel(cl)
# 
# print(paste("Getting plots for", length(my.genes), "genes.", Sys.time()))
# m.list <- foreach(i = 1:length(my.genes),
#                   .packages="ggplot2") %dopar% {
#                     jgene = my.genes[i]
#                     # get number of tissues found rhythmic for jgene: use it in title
#                     count <- r.genes.sum[r.genes.sum$gene == jgene, "count"]
#                     jtitle <- paste(jgene, "rhythmic in", count, "tissues")
#                     m <- PlotAcrossTissues(dat.sub.temp, fit.list, jgene, jtitle)
#                     m
#                   }
# print(paste("Done getting plots.", Sys.time()))
# 
# fit.list <- foreach(i = 1:length(tissues), 
#                     .packages="plyr") %dopar% {
#                       jtissue <- tissues[i]
#                       dat.tiss <- subset(dat, tissue %in% c(jtissue))
#                       fit.out <- dlply(dat.tiss, .(gene), GetAmpPhase)
#                     }
# 
# # print to pdf file
# pdf("plots/rhythmic_genes_across_tissues2.pdf")
# lapply(m.list, print)
# dev.off()
# rm(m.list)