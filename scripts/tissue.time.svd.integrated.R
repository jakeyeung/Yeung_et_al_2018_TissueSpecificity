# tissue.time.svd.integrated.R
# Jake Yeung
# Dec 5 2015
# Combine RNA-Seq wih microarray: do Fourier analysis


# Functions ---------------------------------------------------------------

scripts.dir <- "scripts"
funcs.dir <- file.path(scripts.dir, "functions")
source(file.path(funcs.dir, "DataHandlingFunctions.R"))
source(file.path(funcs.dir, "GetTissueTimes.R"))
source(file.path(funcs.dir, "RemoveProblemGenes.R"))

# Define dirs -------------------------------------------------------------

# define dirs
data.dir <- "data"
unadj.array.fname <- "array_exprs_colnames_fixed.best.probe.selected.txt"
unadj.array.path <- file.path(data.dir, unadj.array.fname)
normalized.array.fname <- "array.adj.0.07.txt"
normalized.array.path <- file.path(data.dir, normalized.array.fname)
rna.seq.fname <- "rna_seq_deseq_counts_colnames_fixed.txt"
rna.seq.path <- file.path(data.dir, rna.seq.fname)


# Load file ---------------------------------------------------------------

unadj.array <- read.table(unadj.array.path, header=TRUE, sep='\t')
normalized.array <- read.table(normalized.array.path)
rna.seq.exprs <- read.table(rna.seq.path, header=TRUE, sep='\t')



# Handle array unadjusted -------------------------------------------------

# first column contained gene names
rownames(unadj.array) <- make.names(unadj.array$gene, unique=TRUE)

drop.cols <- c("gene")
unadj.array <- unadj.array[, !(names(unadj.array) %in% drop.cols)]
Peek(unadj.array)


# Transform to normal scale -----------------------------------------------

unadj.array <- as.data.frame((2^unadj.array))

# Handle array ------------------------------------------------------------

normalized.array <- RemoveProblemGenes(normalized.array)

# Handle duplicate rownames: RNASEQ ------------------------------------

rownames(rna.seq.exprs) <- make.names(rna.seq.exprs$gene, unique=TRUE)

drop.cols <- c("gene")
rna.seq.exprs <- rna.seq.exprs[, !(names(rna.seq.exprs) %in% drop.cols)]

Peek(rna.seq.exprs)  # expect gene names as row names, tissues in columns


# Define constants --------------------------------------------------------

tissues <- GetTissues(colnames(normalized.array)) 

# Define common genes -----------------------------------------------------

filtered.genes <- intersect(rownames(normalized.array), rownames(rna.seq.exprs))

normalized.array <- normalized.array[filtered.genes, ]
rna.seq.exprs <- rna.seq.exprs[filtered.genes, ]
unadj.array <- unadj.array[filtered.genes, ]

# Combine RNA-Seq with Array ----------------------------------------------

# Merge data into long format ---------------------------------------------

tissues <- GetTissues(colnames(normalized.array))
times.array <- GetTimes(colnames(normalized.array))
times.rnaseq <- GetTimes(colnames(rna.seq.exprs))

long.array <- data.frame(gene=rep(filtered.genes, length(tissues) * length(times.array)),
                         tissue=as.factor(rep(tissues, each=(length(filtered.genes) * length(times.array)))),
                         time=as.numeric(rep(times.array, each=length(filtered.genes))),
                         experiment=as.factor("array"),
                         exprs=unlist(normalized.array))

long.array.unadj <- data.frame(gene=rep(filtered.genes, length(tissues) * length(times.array)),
                               tissue=as.factor(rep(tissues, each=(length(filtered.genes) * length(times.array)))),
                               time=as.numeric(rep(times.array, each=length(filtered.genes))),
                               experiment=as.factor("unadj.array"),
                               exprs=unlist(unadj.array))

long.rnaseq <- data.frame(gene=rep(filtered.genes, length(tissues) * length(times.rnaseq)),
                          tissue=rep(tissues, each=(length(filtered.genes) * length(times.rnaseq))),
                          time=as.numeric(rep(times.rnaseq, each=length(filtered.genes))),
                          experiment=as.factor("rnaseq"),
                          exprs=unlist(rna.seq.exprs))

dat <- rbind(long.array, long.rnaseq, long.array.unadj)

str(dat)


# Collect garbage ---------------------------------------------------------

rm(long.array, long.array.unadj, long.rnaseq)


# Plot expressions --------------------------------------------------------


# Load genes from file ----------------------------------------------------

# modules <- seq(5)
# N <- 20  # get top 100 genes
# 
# for (module in modules){
#   fname <- file.path("results", paste0("module_", module, "_genes.txt"))
#   outfile <- file.path("plots", paste0("module_", module, "_genes.exprs.pdf"))
#   pdf(outfile)
#   gene.list <- system(paste("cut -f1", fname), intern = TRUE)[1:N]
#   for (gene in gene.list){
#     dat.sub <- subset(dat, (gene == gene & experiment == "array"))
#     m <- ggplot(dat.sub, aes(x=time, y=log2(exprs + 1), colour=tissue, shape=experiment)) + 
#       geom_point(size=4) + 
#       geom_line() + 
#       ggtitle(paste(mygene)) 
#     print(m)
#   }
#   dev.off()
# }
# 
# # Plot expression of a gene
# # cool genes from module 1 of SVD analysis...
# gene.list <- c("Serpine1", "Lonrf1", "Asb12", "Spon2", "Pnpla3", "Dhrs9", "Cys1", "Tspan4", "Svs5", "Fam124b", "Lcn9", "Rorc", "Crisp1")
# # tissues.sub <- c("Adr", "Aorta", "Kidney", "WFAT")
# tissues.sub <- tissues
# experiments <- c("array")  # array, rnaseq, unadj.array
# 
# for (mygene in gene.list){
#   dat.sub <- subset(dat, (gene==mygene & tissue %in% tissues.sub & experiment %in% experiments))
#   m <- ggplot(dat.sub, aes(x=time, y=log2(exprs + 1), colour=tissue, shape=experiment)) + 
#       geom_point(size=4) + 
#       geom_line() + 
#       ggtitle(paste(mygene)) 
#   print(m)
# }
# 
# # plot top genes from module 2 of SVD analysis
# gene.list <- c("Lcn8", "Gpx5", "Ly6g5b", "Lcn9", "Cst8", "Defb47", "Ly6g5c")
# tissues.sub <- tissues
# experiments <- c("array", "rnaseq")
# 
# for (mygene in gene.list){
#   dat.sub <- subset(dat, (gene==mygene & tissue %in% tissues.sub & experiment %in% experiments))
#   m <- ggplot(dat.sub, aes(x=time, y=log2(exprs + 1), colour=tissue, shape=experiment)) + 
#     geom_point(size=4) + 
#     geom_line() + 
#     ggtitle(paste(mygene)) 
#   print(m)
# }



