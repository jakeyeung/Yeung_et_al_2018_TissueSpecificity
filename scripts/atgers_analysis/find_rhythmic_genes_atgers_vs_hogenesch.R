# 2016-02-21
# Jake Yeung
# Find rhythmic genes in Atger et al

library(dplyr)
library(parallel)
library(ggplot2)
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/PlotFunctions.R")
source("scripts/functions/NcondsFunctions.R")
source("scripts/functions/Queue.R")
source("scripts/functions/ListFunctions.R")
source("scripts/functions/MemoryManagement.R")
source("scripts/functions/EnvironmentsFunctions.R")

inf <- "/home/yeung/data/tissue_specificity/atger_et_al/GSE73554_WT_AL_Intron_Exon_RFP.txt"
dat <- read.delim2(inf, header = TRUE, sep = "\t")

load("Robjs/dat.long.fixed_rik_genes.Robj")




# Keep columns with exons -------------------------------------------------

keep.cols <- "*Exon*"

dat.filt <- dat[, grepl(keep.cols, colnames(dat)), ]

rownames(dat.filt) <- make.unique(as.character(dat$Gene_Symbol))

times <- as.numeric(sapply(colnames(dat.filt), function(cname) strsplit(cname, "_")[[1]][[4]]))
samps <- sapply(colnames(dat.filt), function(cname) strsplit(cname, "_")[[1]][[5]])

dat.atger.long <- data.frame(gene = dat$Gene_Symbol,
                             tissue = "Liver",
                             time = rep(times, each = nrow(dat.filt)),
                             experiment = "atger",
                             exprs = as.numeric(as.character(unlist(dat.filt))),
                             samp = rep(samps, each = nrow(dat.filt)))

dat.atger.long$samp.numeric <- as.numeric(chartr(old = "ABCD", new = '1234', x = dat.atger.long$samp))

PlotGeneAcrossTissues(subset(dat.atger.long, gene == "Gnai3"))


# Stagger samples across 96 hours  ----------------------------------------

dat.atger.long$time.split <- apply(dat.atger.long, 1, function(row){
  time <- as.numeric(row[3]); samp.num <- as.numeric(row[7])
  return(time + 24 * (samp.num - 1))
})

dat.atger.long$time <- dat.atger.long$time.split
dat.atger.long$samp <- NULL
dat.atger.long$samp.numeric <- NULL
dat.atger.long$time.split <- NULL

common.genes <- intersect(unique(as.character(dat.long$gene)), unique(as.character(dat.atger.long$gene)))
dat.merged <- rbind(subset(dat.long, tissue == "Liver" & gene %in% common.genes), subset(dat.atger.long, gene %in% common.genes))

PlotGeneAcrossTissues(subset(dat.atger.long, gene == "Tars"))

PlotGeneAcrossTissues(subset(dat.merged, gene == "Dbp"))


# Setup data for Nconds ---------------------------------------------------

dat.merged$tissue <- sapply(dat.merged$experiment, function(ex){
  if (ex == "array" | ex == "rnaseq"){
    return("LiverHogenesch")
  } else if (ex == "atger"){
    return("LiverAtger")
  }
})

PlotGeneAcrossTissues(subset(dat.merged, gene == "Dbp"))

print("Chunking data to environment")
dat.env <- DatLongToEnvironment(dat.merged)

start <- Sys.time()
# Rprof()
fits.all <- mclapply(ls(dat.env), function(gene){
  MakeDesMatRunFitEnv(dat.env, gene, unique(dat.merged$tissue), n.rhyth.max = 2, w = 2 * pi / 24, criterion = "BIC", normalize.weights = TRUE, cutoff = 1e-5, top.n = NULL, sparse = FALSE)
}, mc.cores = 5)
print(Sys.time() - start)

fits.all.long <- lapply(fits.all, function(x){
  gene <- x$gene
  x$gene <- NULL
  fits.temp.long <- ListToLong(x, gene, top.n = 5)
})
fits.all.long <- do.call(rbind, fits.all.long)
print(Sys.time() - start)

fits.all.long$n.params <- sapply(fits.all.long$model, function(m) return(length(strsplit(as.character(m), ";")[[1]])))
fits.all.long$n.rhyth <- sapply(fits.all.long$model, GetNrhythFromModel)
fits.all.long$amp.avg <- mapply(GetAvgAmpFromParams, fits.all.long$param.list, fits.all.long$model)
fits.all.long$phase.sd <- mapply(GetSdPhaseFromParams, fits.all.long$param.list, fits.all.long$model)
fits.all.long$phase.maxdiff <- mapply(GetMaxPhaseDiffFromParams, fits.all.long$param.list, fits.all.long$model)

# get best model from each n.params
fits.best.nparam <- fits.all.long %>%
  group_by(gene) %>%
  filter(weight.raw == min(weight.raw))

save(fits.all.long, file = "Robjs/atger_vs_hogenesch.fits.all.long.Robj")
save(fits.best.nparam, file = "Robjs/atger_vs_hogenesch.fits.all.long.Robj")

# How many in each model? -------------------------------------------------

n.models <- fits.best.nparam %>%
  group_by(model) %>%
  summarise(n.models = length(model))

n.models$model <- as.character(n.models$model)
n.models$model <- as.factor(c("Flat", "Hogenesch", "Atger", "Hog-At-Shared", "Hog-At-Different"))
n.models <- OrderDecreasing(n.models, jfactor = "model", jval = "n.models")
ggplot(n.models, aes(x = model, y = n.models)) + geom_bar(stat = "identity") + theme_bw(24) + 
  ggtitle("Nconds: Atger (LD) vs Hogenesch (DD) Liver rhythmic models") + xlab("") + ylab("Number of models")

# many genes in LiverHogenesch;LiverAtger model, are the amplitudes large?

jtest <- subset(fits.best.nparam, model == "LiverHogenesch;LiverAtger")
head(jtest[order(jtest$amp.avg, decreasing = TRUE), ])
PlotGeneAcrossTissues(subset(dat.merged, gene == "Mfsd2a"))
PlotGeneAcrossTissues(subset(dat.merged, gene == "Npas2"))
PlotGeneAcrossTissues(subset(dat.merged, gene == "Dhrs9"))
PlotGeneAcrossTissues(subset(dat.merged, gene == "Wisp2"))
PlotGeneAcrossTissues(subset(dat.merged, gene == "Car5a"))
PlotGeneAcrossTissues(subset(dat.merged, gene == "Elovl3"))


jtest <- subset(fits.best.nparam, model == "LiverAtger")
head(jtest[order(jtest$amp.avg, decreasing = TRUE), ])
PlotGeneAcrossTissues(subset(dat.merged, gene == "Moxd1"))
PlotGeneAcrossTissues(subset(dat.merged, gene == "Tnnc1"))
PlotGeneAcrossTissues(subset(dat.merged, gene == "Lmod2"))
PlotGeneAcrossTissues(subset(dat.merged, gene == "Camta1"))

jtest <- subset(fits.best.nparam, model == "LiverHogenesch")
head(jtest[order(jtest$amp.avg, decreasing = TRUE), ])
PlotGeneAcrossTissues(subset(dat.merged, gene == "Cyp26a1"))
PlotGeneAcrossTissues(subset(dat.merged, gene == "Tst"))
