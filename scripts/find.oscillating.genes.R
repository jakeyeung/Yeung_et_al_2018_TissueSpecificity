# Jake Yeung
# find oscillating genes
# find.oscillating.genes.R
# Dec 4 2014

# Adjustable parameters ---------------------------------------------------

outfile <- "data/oscillating_genes.txt"
gene.header <- "gene"
params.header <- c("amp", "phase", "int.array", "int.rnaseq")

# Define constants --------------------------------------------------------

pval.threshold <- 1e-5
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
FitComplexModel <- function(df){
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
  amp <- sqrt(params[["a"]] ^ 2 + params[["b"]] ^ 2)
  phase <- atan(params[["b"]] / params[["a"]])
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

fit.list <- vector(mode="list", length=length(tissues))
names(fit.list) <- tissues


# Takes ~30 minutes
for (jtissue in tissues){
  print(Sys.time())
  print(paste0("Fitting models for tissue: ", jtissue))
  # dat.tiss <- subset(dat, tissue %in% c(jtissue) & gene %in% clockgenes)
  dat.tiss <- subset(dat, tissue %in% c(jtissue))
  # fit.out <- dlply(dat.tiss, .(gene), FindRhythmic)
  fit.out <- dlply(dat.tiss, .(gene), GetAmpPhase)
  fit.list[[jtissue]] <- fit.out
}


# Write list to textfile --------------------------------------------------

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


# Plot distribution of amplitudes -----------------------------------------

GetFitInfo(fit.list[[tissue]])

amplitudes <- vector(mode="list", length(tissues))
names(amplitudes) <- tissues

for (tissue in tissues){
  amplitudes.temp <- lapply(fit.list[[tissue]], function(x){
    return(x[["amp"]])
  })
  amplitudes[[tissue]] <- amplitudes.temp
}
rm(amplitudes.temp)


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

ggplot(subset(fits.df, pval <= 1e-5), aes(x = amp, fill = tissue)) + 
  geom_density() + 
  facet_wrap(~tissue) + 
  scale_x_log10()


# Get top rhythmic genes --------------------------------------------------

fits.df.subset <- subset(fits.df, pval < 1e-5)

# order by tissue, then by decreasing amp
fits.df.subset <- fits.df.subset[order(fits.df.subset$tissue, 
                                       -fits.df.subset$amp), ]

# Get top 20 genes (by amplitude) for each tissue
top.n <- 20

fits.df.subset.top.n <- data.frame(tissue = character(),
                                   gene = character(),
                                   amp = character(),
                                   phase = character(),
                                   pval = character())

for (jtissue in tissues){
  tissue.temp <- head(subset(fits.df.subset, tissue == jtissue), n=top.n)
  fits.df.subset.top.n <- rbind(fits.df.subset.top.n, tissue.temp)
}

head(fits.df.subset.top.n)
tail(fits.df.subset.top.n)

# summmarise
ddply(fits.df.subset.top.n, .(tissue), summarise, 
      count = length(amp), 
      avg = mean(amp), 
      min = min(amp),
      max = max(amp))


# Get top oscillating genes -----------------------------------------------

# Do it by tissues

rhythmic.genes <- vector(mode="list", length(tissues))
names(rhythmic.genes) <- tissues

for (tissue in tissues){
  # stolen from: http://stackoverflow.com/questions/14771029/filter-values-from-list-in-r
  rhythmic.genes.temp <- lapply(fit.list[[tissue]], GetRhythmicGene, 
                                pval.threshold, amp.threshold=0.1) == TRUE
  # Keep only element that are "TRUE"
  rhythmic.genes.temp <- names(fit.list[[tissue]])[rhythmic.genes.temp]
  rhythmic.genes[[tissue]] <- rhythmic.genes.temp
}


# THIS IS SUPER SLOW.
# for (gene in c(clockgenes, tissuegenes)){
#   dat.gene.tiss <- subset(dat[dat$gene == gene, ], tissue=="Adr")
#   dat.gene.tiss <- dat.gene.tiss[order(dat.gene.tiss$time), ]
#   # all negatives go to 0
#   dat.gene.tiss$exprs[which(dat.gene.tiss$exprs < 0)] <- 0
#   fit.complex <- lm(exprs ~ 0 + experiment + cos(omega * time) + sin(omega * time), 
#                     data = dat.gene.tiss)
#   # low pvalue, it is worth it to use a complex model
#   fit <- fit.complex
#   # get pval for complex model by comparing with model containing only intercepts
#   fit.flat <- lm(exprs ~ 0 + experiment, data = dat.gene.tiss)
#   ftest.flat <- anova(fit.flat, fit.complex)
#   pval <- ftest.flat["Pr(>F)"][[1]][2]
#   fit.list[[gene]] <- list(fit=fit, pval=pval)
# }


# Get top genes -----------------------------------------------------------

pval.df <- data.frame(pvals=rep(NA, length(filtered.genes)), gene=filtered.genes, row.names=filtered.genes)

for (gene in c(clockgenes, tissuegenes)){
  pval <- fit.list[[gene]][["pval"]]
  pval.df[gene, "pvals"] <- pval 
}

pval.df <- na.omit(pval.df)

# Plot top genes ----------------------------------------------------------

pval.df.ordered <- pval.df[with(pval.df, order(pvals)), ]
top.genes <- pval.df.ordered$gene

pdf(fit.plot)
for (gene in top.genes){
  pval <- pval.df.ordered[gene, ]$pvals
  dat.gene.tiss <- subset(dat[dat$gene == gene, ], tissue=="Adr")
  dat.gene.tiss <- dat.gene.tiss[order(dat.gene.tiss$time), ]
  
  fit <- fit.list[[gene]][["fit"]]
  
  # Get vectors for plotting
  y.combined <- dat.gene.tiss$exprs
  t.combined <- dat.gene.tiss$time
  y.array <- dat.gene.tiss[dat.gene.tiss$experiment=="array", ]$exprs
  t.array <- dat.gene.tiss[dat.gene.tiss$experiment=="array", ]$time
  y.rnaseq <- dat.gene.tiss[dat.gene.tiss$experiment=="rnaseq", ]$exprs
  t.rnaseq <- dat.gene.tiss[dat.gene.tiss$experiment=="rnaseq", ]$time
  
  # Plot GIVENS
  plot(t.combined, y.combined, main=paste(gene, "Pval:", signif(pval, 2)))
  points(t.rnaseq, y.rnaseq, pch='*')
  
  # BEGIN: get and plot y.models
  t <- seq(18, 64, length.out = 100)  # tspan
  # get fit parameters: first ask if simple or complex
  if (IsSimple(fit)){
    y.model.combined <- oscillate(params = coef(fit), omega, t)  
    lines(t, y.model.combined, col='red')
  } else {
    # take 1st intercept, 3 and 4 are cos and sin coefficients
    y.model.array <- oscillate(params = coef(fit)[c(1, 3, 4)], omega, t)
    y.model.rnaseq <- oscillate(params = coef(fit)[c(2, 3, 4)], omega, t)    
    lines(t, y.model.array, col='blue')
    lines(t, y.model.rnaseq, col='orange')
  }
  legend("topleft", c("combined.fit", "array.fit", "rnaseq.fit"),
         lty=c(1,1),
         lwd=c(2.5, 2.5),
         col=c("red", "blue", "orange"))  
}
dev.off()
