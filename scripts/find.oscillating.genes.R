# Jake Yeung
# find oscillating genes
# find.oscillating.genes.R
# Dec 4 2014


# Define constants --------------------------------------------------------

pval.threshold <- 1e-3
T <- 24  # hours
omega <- 2 * pi / T


# Functions ---------------------------------------------------------------

scripts.dir <- 'scripts'
funcs.dir <- 'functions'
source(file.path(scripts.dir, funcs.dir, "ReplaceNegs.R"))
source(file.path(scripts.dir, funcs.dir, "GetTissueTimes.R"))  # get tissue and times
source(file.path(scripts.dir, funcs.dir, "DataHandlingFunctions.R"))
source(file.path(scripts.dir, funcs.dir, "RegressionFunctions.R"))

# oscillating function
oscillate <- function(params, omega, t){
  params[1] + params[2] * cos(omega * t) + params[3] * sin(omega * t)
}


# Define dirs -------------------------------------------------------------

# define dirs
data.dir <- "data"
normalized.array.fname <- "array.adj.0.07.txt"
normalized.array.path <- file.path(data.dir, normalized.array.fname)
rna.seq.fname <- "rna_seq_deseq_counts_colnames_fixed.txt"
rna.seq.path <- file.path(data.dir, rna.seq.fname)

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

clockgenes <- c("Nr1d1","Dbp","Arntl","Npas2","Nr1d2","Bhlhe41","Nfil3",
                "Cdkn1a","Lonrf3","Tef","Usp2","Wee1","Dtx4","Asb12",
                "Elovl3","Clock","Per1","Per2","Per3","Cry2","Cry1")
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

fit.list <- vector(mode="list", length=length(filtered.genes))
names(fit.list) <- filtered.genes

for (gene in clockgenes){
  dat.gene.tiss <- subset(dat[dat$gene == gene, ], tissue=="Adr")
  dat.gene.tiss <- dat.gene.tiss[order(dat.gene.tiss$time), ]
  # all negatives go to 0
  dat.gene.tiss$exprs[which(dat.gene.tiss$exprs < 0)] <- 0
  fit.complex <- lm(exprs ~ 0 + experiment + cos(omega * time) + sin(omega * time), 
                    data = dat.gene.tiss)
  fit.simple <- lm(exprs ~ cos(omega * time) + sin(omega * time), 
                   data = dat.gene.tiss)
  
  # model selection
  ftest <- anova(fit.simple, fit.complex)
  ftest.pval <- ftest["Pr(>F)"][[1]][2]
  
  if (ftest.pval < pval.threshold) {
    # low pvalue, it is worth it to use a complex model
    fit <- fit.complex
  } else
    # stick with simple model
    fit <- fit.simple
  
  pval <- PvalFromFit(fit)
  
  fit.list[[gene]] <- list(fit=fit, pval=pval)
  
}


# Get top genes -----------------------------------------------------------


# Plot top genes ----------------------------------------------------------

# top.genes <- 
  
for (gene in top.genes){
  dat.gene.tiss <- subset(dat[dat$gene == gene, ], tissue=="Adr")
  dat.gene.tiss <- dat.gene.tiss[order(dat.gene.tiss$time), ]
  # Get vectors for plotting
  y.combined <- dat.gene.tiss$exprs
  t.combined <- dat.gene.tiss$time
  y.array <- dat.gene.tiss[dat.gene.tiss$experiment=="array", ]$exprs
  t.array <- dat.gene.tiss[dat.gene.tiss$experiment=="array", ]$time
  y.rnaseq <- dat.gene.tiss[dat.gene.tiss$experiment=="rnaseq", ]$exprs
  t.rnaseq <- dat.gene.tiss[dat.gene.tiss$experiment=="rnaseq", ]$time
  
  t <- seq(18, 64, length.out = 100)
  y.model.combined <- oscillate(params = coef(fit.simple), omega, t) 
  # take 1st intercept, 3 and 4 are cos and sin coefficients
  y.model.array <- oscillate(params = coef(fit.complex)[c(1, 3, 4)], omega, t)
  y.model.rnaseq <- oscillate(params = coef(fit.complex)[c(2, 3, 4)], omega, t)
  
  # Plot vectors
  plot(t.combined, y.combined, main=gene)
  points(t.rnaseq, y.rnaseq, pch='*')
  lines(t, y.model.combined, col='red')
  lines(t, y.model.array, col='blue')
  lines(t, y.model.rnaseq, col='orange')
  legend("topleft", c("combined.fit", "array.fit", "rnaseq.fit"),
         lty=c(1,1),
         lwd=c(2.5, 2.5),
         col=c("red", "blue", "orange"))  
}

                                   