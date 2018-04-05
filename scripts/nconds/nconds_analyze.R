# Jake Yeung
# 2015-09-24
# Split ncond outputs into "tissue-specific" and "tissue-wide"

library(ggplot2)
library(hash)

NumberRhythmicModel <- function(matobj, model){
  # get number of rhythmic given model number
  cnames <- colnames(matobj[[model]])
  cnames <- cnames[grepl("a.", cnames)]
  n_tissues <- length(strsplit(paste0(cnames[grepl("a", cnames)], collapse = ","), ",")[[1]])
  return(n_tissues)
}

NumberRhythmic <- function(fitobj, matobj, gene){
  cnames <- colnames(matobj[[fitobj[[gene]][["model"]]]])
  # colnames containing a. or b. suggests it is rhythmic in tissue
  cnames <- cnames[grepl("a.", cnames)]
  # get length number of tissues that are rhythmic by counting how many a.#s there are
#   print(paste0(cnames[grepl("a.", cnames)], collapse = ","))
#   print(strsplit(paste0(cnames[grepl("a.", cnames)], collapse = ","), ",")[[1]])
  n_tissues <- length(strsplit(paste0(cnames[grepl("a", cnames)], collapse = ","), ",")[[1]])
  return(n_tissues)
}

WhichRhythmic <- function(fitobj, matobj, gene){
  model <- fitobj[[gene]][["model"]]
  cnames <- colnames(matobj[[model]])
  cnames <- cnames[grepl("a.", cnames)]
  tissuesi <- strsplit(paste0(cnames[grepl("a", cnames)], collapse = ","), ",")[[1]]
  tissues.i <- sapply(tissuesi, function(s) as.numeric(strsplit(s, "a.")[[1]][[2]]), USE.NAMES = FALSE)
  if (length(tissues.i) == 0){
    return("")
  }
  tissues <- unique(rownames(matobj[[model]]))
  rhyth.tissues <- tissues[tissues.i]
  return(rhyth.tissues)
}

# Load --------------------------------------------------------------------

outobj <- "/home/yeung/projects/tissue-specificity/nconds_outputs/nconds_8_conds_rerun_less_cores/nconds_8_tissues.fit_output.Robj"
outmat <- "/home/yeung/projects/tissue-specificity/nconds_outputs/nconds_8_conds_rerun_less_cores/nconds_8_tissues.matrix_output.Robj"
load(outobj, verbose = T)
load(outmat, verbose = T)

load("Robjs/dat.fit.Robj", verbose = T)

# Find models that are "tissue-specific" versus "tissue-wide" -------------

WhichRhythmic(fit, my_mat, "Xrcc1")

n_rhythmic <- lapply(names(fit), function(gene){
  n <- NumberRhythmic(fit, my_mat, gene)
  return(data.frame(n.rhyth = n, gene = gene, model = fit[[gene]][["model"]]))
})
dat.nrhyth <- do.call(rbind, n_rhythmic)

dat.nrhyth$rhyth.tiss <- sapply(dat.nrhyth$gene, function(gene){
  rhyth.tiss <- WhichRhythmic(fit, my_mat, gene)
  return(paste0(rhyth.tiss, collapse = ","))
})


# Get average and variance of amplitude of rhythmic gene ------------------

tissue.amp.dic <- hash(paste(dat.fit$gene, dat.fit$tissue, sep = ","), dat.fit$amp)

GetRhythmicProperties <- function(gene, rhyth.tiss, tissue.amp.dic, statistic = "mean"){
  rhyth.tiss.vec <- strsplit(rhyth.tiss, split = ",")[[1]]
  vec <- numeric(length = length(rhyth.tiss.vec))
  for (i in 1:length(vec)){
    tiss <- rhyth.tiss.vec[i]
    key <- paste(gene, tiss, sep = ",")
    amp <- tissue.amp.dic[[key]]
    if (is.null(amp)){
      return(NA)
    }
    vec[i] <- tissue.amp.dic[[key]]
  }
  if (statistic == "mean"){
    return(mean(vec))
  } else if (statistic == "variance"){
    return(var(vec))
  }
}

dat.nrhyth$amp.var <- mapply(GetRhythmicProperties, 
                             dat.nrhyth$gene, dat.nrhyth$rhyth.tiss,
                             MoreArgs=list(tissue.amp.dic = tissue.amp.dic, statistic = "variance"))

dat.nrhyth$amp.avg <- mapply(function(gene, rhyth.tiss){
#   gene <- xrow[2]
#   rhyth.tiss <- xrow[4]
  rhyth.tiss.vec <- strsplit(rhyth.tiss, split = ",")[[1]]
  vec <- numeric(length = length(rhyth.tiss.vec))
  for (i in 1:length(vec)){
    tiss <- rhyth.tiss.vec[i]
    key <- paste(gene, tiss, sep = ",")
    amp <- tissue.amp.dic[[key]]
    if (is.null(amp)){
      return(NA)
    }
    vec[i] <- tissue.amp.dic[[key]]
  }
  return(mean(vec))
  }, dat.nrhyth$gene, dat.nrhyth$rhyth.tiss)

# Plot distribution -------------------------------------------------------

ggplot(subset(dat.nrhyth, n.rhyth > 0), aes(x = n.rhyth)) + geom_histogram()

dat.nrhyth.sub <- subset(dat.nrhyth, n.rhyth == 4)
dat.nrhyth.sub[order(dat.nrhyth.sub$model), ]

dat.n.rhyth.subcount <- dat.nrhyth.sub %>% group_by(model) %>% summarise(count = length(model))
dat.n.rhyth.subcount[order(dat.n.rhyth.subcount$count, decreasing = TRUE), ]  # model 571 is most

subset(dat.nrhyth, n.rhyth == 4 & model == 571)

jgene <- "Znhit6"
jgene <- "Zhx2"  # rhythmic
jgene <- "Plek"
jgene <- "Klf6"
jgene <- "Usp2"
jgene <- "Zmynd11"
jgene <- "Fam69b"
jgene <- "Fh1"
jgene <- "Rps23"
WhichRhythmic(fit, my_mat, jgene)
PlotGeneAcrossTissues(subset(dat.long, gene == jgene))

# save(dat.nrhyth, file = "Robjs/nconds.dat.nrhyth.Robj")


# Number of models with number of tissues ---------------------------------

models <- seq(length(my_mat))

n.rhyth.bg <- lapply(models, function(i){
  n <- NumberRhythmicModel(my_mat, i)
  return(data.frame(model.i = i, n = n))
})

n.rhyth.bg.df <- do.call(rbind, n.rhyth.bg)

n.rhyth.bg.sum <- n.rhyth.bg.df %>%
  group_by(n) %>%
  summarise(count = length(n))

ggplot(n.rhyth.bg.sum, aes(x = n, y = count)) + geom_bar(stat = "identity")
