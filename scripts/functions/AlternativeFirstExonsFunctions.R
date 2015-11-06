DoCanCor <- function(tpm.sub, xvar = "tpm_norm.avg", yvar = "amp"){
  # do canonical correlation of two matrices to find interesting
  # alternative promoter usage
  if (nrow(tpm.sub) == 0){
    return(data.frame(Cor = NA))
  }
  X <- dcast(data = tpm.sub, formula = transcript_id ~ tissue, value.var = xvar)
  rownames(X) <- X$transcript_id; X$transcript_id <- NULL
  jtrans <- tpm.sub$transcript_id[[1]]
  Y <- data.frame(subset(tpm.sub, transcript_id == jtrans, select = c("amp", "tissue")))
  rownames(Y) <- as.character(Y$tissue); Y$tissue <- NULL
  jcor <- cancor(tpm.afe.mat, dat.fit.mat, xcenter = F, ycenter = F)
  return(data.frame(Cor = jcor$cor))
}

KeepUpToThres <- function(vec, thres, min.dim = 2){
  # keep vector up to threshold return by index
  for (i in seq(length(vec))){
    if (vec[i] > thres){
      break
    }
  }
  if (i < min.dim){
    return(i:min.dim)
  }
  return(1:i)
}

GetPromoterUsage <- function(dat, do.svd = TRUE, thres = 0.9, append.tiss = TRUE, get.means=TRUE, get.entropy=TRUE){
  # get promoter usage
  dat.mat <- dcast(dat, tissue + amp + mean~ transcript_id, value.var = "tpm_norm.avg")
  dat.mat.prom <- subset(dat.mat, select = -c(tissue, amp, mean))
  
  dat.H <- apply(dat.mat.prom, 1, ShannonEntropy)
  if (do.svd){
    # reduce dim
    dat.mat.prom <- sweep(dat.mat.prom, MARGIN = 1, STATS = rowMeans(dat.mat.prom), FUN = "-")
    dat.mat.prom.s <- svd(dat.mat.prom)
    eigvals <- dat.mat.prom.s$d ^ 2 / sum(dat.mat.prom.s$d ^ 2)
    eigvals.cum <- cumsum(eigvals)
    # threshold for number of eigvals to keep
    keep <- KeepUpToThres(eigvals.cum, thres, min.dim = 2)
    dat.mat.prom.trans <- sweep(dat.mat.prom.s$u[, keep], MARGIN = 2, STATS = dat.mat.prom.s$d[keep], FUN = "*")
    dat.mat.trans <- data.frame(amp = dat.mat$amp, dat.mat.prom.trans) 
  } else {
    dat.mat.trans <- data.frame(amp = dat.mat$amp, dat.mat.prom)
  }
  if (append.tiss){
    dat.mat.trans$tissue <- dat.mat$tissue
  }
  out <- list(dat.mat.trans = dat.mat.trans)
  if (get.means){
    out$weights <- dat.mat$mean
  }
  if (get.entropy){
    out$entropy <- dat.H
  }
  return(out)
}

CorrelateAmpPromMulti <- function(dat, thres = 0.9, do.svd = TRUE, weighted = FALSE, eps = 1e-10){
  #   dat.mat <- dcast(dat, tissue + amp + mean~ transcript_id, value.var = "tpm_norm.avg")
  #   dat.mat.prom <- subset(dat.mat, select = -c(tissue, amp, mean))
  #   
  #   # reduce dim
  #   dat.mat.prom <- sweep(dat.mat.prom, MARGIN = 1, STATS = rowMeans(dat.mat.prom), FUN = "-")
  #   dat.mat.prom.s <- svd(dat.mat.prom)
  #   eigvals <- dat.mat.prom.s$d ^ 2 / sum(dat.mat.prom.s$d ^ 2)
  #   eigvals.cum <- cumsum(eigvals)
  #   # threshold for number of eigvals to keep
  #   keep <- KeepUpToThres(eigvals.cum, thres, min.dim = 2)
  #   dat.mat.prom.trans <- dat.mat.prom.s$u[, keep] * dat.mat.prom.s$d[keep]
  #   
  #   dat.mat.trans <- data.frame(amp = dat.mat$amp, dat.mat.prom.trans)
  dat.mat.trans.lst <- GetPromoterUsage(dat, do.svd = do.svd, thres = thres, append.tiss = FALSE, get.means = TRUE, get.entropy = TRUE)
  dat.mat.trans <- dat.mat.trans.lst$dat.mat.trans
  weights <- dat.mat.trans.lst$weights
  entropy <- dat.mat.trans.lst$entropy
  if (weighted){
    # by shannon entropy
    weights <- 1 - entropy
    weights <- weights + eps  # if zero?
    fit.altprom <- lm(formula = amp ~ ., data = dat.mat.trans, weights = weights)
  } else {
    fit.altprom <- lm(formula = amp ~ ., data = dat.mat.trans)
  }
  # return just R-squared and best p-value
  rsqr <- summary(fit.altprom)$r.squared
  coefs <- summary(fit.altprom)$coefficients
  pval.best <- min(coefs[2:nrow(coefs), "Pr(>|t|)"])
  return(data.frame(rsqr = rsqr, pval.best = pval.best))
}

# Functions: kallisto -----------------------------------------------------

PlotDiagnostics <- function(dat.tpm, dat.arrayrnaseq, jgene, jtranscript){
  source('scripts/functions/PlotGeneAcrossTissues.R')
  # uses data objects from find_alternative_first_exons_kallisto.R OR filter_kallisto_for_start_sites.R
  # Plot diagnostic plots, given a gene and transcript
  # jgene <- "Sorbs1"; jtranscript="ENSMUST00000165469"  # Liver and BFAT looks rhythmic, but assigned as "not rhythmic"
  
  test.gene <- subset(dat.tpm, gene_name == jgene)
  test <- subset(test.gene, transcript_id == jtranscript)
  
  rhythmic.tissues <- paste(unique(test$tissue[which(test$is.rhythmic == TRUE)]), collapse = ",")
  notrhythmic.tissues <- paste(unique(test$tissue[which(test$is.rhythmic == FALSE)]), collapse = ",")
  
  print(PlotGeneAcrossTissues(subset(dat.arrayrnaseq, gene == jgene), 
                              jtitle = paste(jgene, "\nRhythmic Tissues:", rhythmic.tissues, "\nFlat Tissues:", notrhythmic.tissues)))
  
  print(PlotTpmAcrossTissues(test.gene, jtitle = paste(jgene, jtranscript), log2.transform=FALSE))
  print(PlotTpmAcrossTissues(test, jtitle = paste(jgene, jtranscript), log2.transform=TRUE))
  
  print(ggplot(test.gene, aes(y = tpm_normalized, x = tissue)) + 
          geom_boxplot() + 
          facet_wrap(~transcript_id) + 
          ggtitle(jgene) +
          theme(axis.text.x = element_text(angle = 90, hjust = 1)))
  
  print(ggplot(test, 
               aes(y = tpm_normalized, x = is.rhythmic)) + 
          geom_boxplot() + 
          ggtitle(paste(jgene, jtranscript)))
}

PlotDiagnostics2 <- function(dat.tpm, tpm.avg, dat.long, jgene, jtranscript){
  source('scripts/functions/PlotGeneAcrossTissues.R')
  # updated version plotting the linear trend
  dat.gene <- subset(dat.tpm, gene_name == jgene)
  dat.transcript <- subset(dat.gene, transcript_id == jtranscript)
  tpm.avg.sub <- subset(tpm.avg, transcript_id == jtranscript)
  
  print(PlotGeneAcrossTissues(subset(dat.long, gene == jgene), jtitle = paste(jgene, jtranscript, sep = "\n")))
  print(PlotTpmAcrossTissues(dat.gene, jtitle = paste("linear-scale", jgene, jtranscript, sep = "\n"), log2.transform=FALSE))
  print(PlotTpmAcrossTissues(dat.transcript, jtitle = paste("log-scale", jgene, jtranscript, sep = "\n"), log2.transform=TRUE))
  print(ggplot(tpm.avg.sub, aes(y = tpm_norm.avg, x = relamp, label = tissue)) + geom_point() + geom_text() + ggtitle(jgene) + geom_smooth(method = "lm"))
}

IsTissueSpecific2 <- function(jdf, pval.min = 1e-5, pval.max = 0.05, cutoff = 4.89){
  # given list of pvals, check if it contains pvals less than pval.min and greater than pval.max.
  # check if pval contains values less than pval.min (rhythmic) and greater than pval.max (flat)
  # if so, return TRUE, otherwise FALSE
  avg.exprs <- jdf$int.rnaseq
  pvals <- jdf$pval
  pvals.filt <- pvals[which(!is.nan(pvals) & avg.exprs > cutoff)]
  if (min(pvals.filt) < pval.min & max(pvals.filt) > pval.max){
    jdf$is.tissue.spec.circ <- rep(TRUE, length(pvals))
  } else {
    jdf$is.tissue.spec.circ<- rep(FALSE, length(pvals))
  }
  return(jdf)
}

IsRhythmic2 <- function(pval, avg.exprs, pval.min = 1e-5, cutoff = 6){
  # ask if gene is rhythmic, not rhythmic or NA (lowly expressed or 0)
  # check if pval is less than pval.min
  if (is.nan(pval) | avg.exprs < cutoff){
    return(NA)
  }
  if (pval < pval.min){
    return(TRUE)
  } else {
    return(FALSE)
  }
}

ModelRhythmicity <- function(dat, jformula=tpm_normalized ~ is.rhythmic){
  # check it contains both Rhythmic and NotRhythmic elements
  if (length(unique(dat$is.rhythmic)) != 2){
    return(data.frame(int = NA, coef = NA, pval = NA))
  }
  jfit <- lm(data = dat, formula = jformula)
  # int <- summary(jfit)$mat["(Intercept)", "Coef"]
  # jcoef <- summary(jfit)$mat["rhythmic.or.notRhythmic", "Coef"]
  int <- coef(jfit)[[1]]
  jcoef <- coef(jfit)[[2]]
  pval <- summary(jfit)$coefficients["is.rhythmicTRUE", "Pr(>|t|)"]
  return(data.frame(int = int, coef = jcoef, pval = pval))
}

CalculateFractionIsoformUsage <- function(tpm, pseudocount = 0){
  # given tpm of a sample across known isoforms, compute
  # fraction of isoform usage. 
  # 
  # Add pseudo count to increase robustness to lowly expressed genes
  tpm_normalized <- (tpm + pseudocount) / (sum(tpm + pseudocount))
}


# Functions: alternative first exons --------------------------------------


KsTestPosNeg <- function(jdf){
  # 
  x.neg <- jdf$sitecount[which(jdf$hit.or.not == "Neg")]
  x.pos <- jdf$sitecount[which(jdf$hit.or.not == "Pos")]
  jks <- ks.test(x.neg, x.pos)
  stat <- jks$statistic
  pval <- jks$p.value
  return(data.frame(statistic = stat, pval = pval))
}

FitPosNeg <- function(jdf){
  jfit <- lm(formula = sitecount ~ hit.or.not, data = jdf)
  summary(jfit)$coefficients["hit.or.notPos", "Pr(>|t|)"]
  pval <- summary(jfit)$coefficients["hit.or.notPos", "Pr(>|t|)"]
  coef <- summary(jfit)$coefficients["hit.or.notPos", "Estimate"]
  return(data.frame(coef = coef, pval = pval))
}

HitOrNot <- function(coef, pval, max.pval=1e-5, min.pval = 0.05){
  if (pval < max.pval){
    if (coef > 0){
      s = "Pos"
    } else if (coef < 0){
      s = "Neg"
    } else {
      warning("Coefficient is 0 but called a hit")
    }
  } else if (pval < min.pval){
    s = NA
  } else {
    s = "NotHit"
  }
  return(s)
}

FitDfToMatrix <- function(jdf, common.genes){
  dat.fitrhyth.filt <- data.frame(subset(jdf, gene %in% common.genes, select = c(gene, tissue, as.numeric(pval), amp)))  # faster
  rnames <- apply(dat.fitrhyth.filt, 1, function(x) paste0(x[2], '-', x[1]))
  rownames(dat.fitrhyth.filt) <- rnames; rm(rnames)
  dat.fitrhyth.filt <- subset(dat.fitrhyth.filt, select = c(pval, amp))
  dat.fitrhyth.filt <- data.matrix(dat.fitrhyth.filt)  # faster to work with matrices
  return(dat.fitrhyth.filt)
}

GetRhythmicOrNot <- function(x, fitdf){
  # Expect x to be a row from cov.normreads, with tissue and gene in 2nd and 3rd col
  tiss <- x[2]
  gene <- x[3]
  rname <- paste0(tiss, '-', gene)
  fitdf.sub = tryCatch({
    fitdf[rname, ]
  }, warning = function(w) {
    print("Warning")
    print(w)
  }, error = function(e) {
    # print(paste("Cannot access:", rname))
    return(NA)
  })  
  if (is.na(fitdf.sub[1])){
    return(NA)
  }
  pval <- fitdf.sub[1]
  amp <- fitdf.sub[2]
  annots <- RhythmicOrNot(pval, amp)
  return(annots)
}

RhythmicOrNot <- function(pval, amp, min.pval = 1e-5, max.pval = 0.05, max.amp = 0.5, min.amp = 0.1){
  if (pval < min.pval & amp > max.amp){
    return("Rhythmic")
  } else if (pval > max.pval & amp < min.amp){
    return("NotRhythmic")
  } else {
    return(NA)  # undecided
  }
}

FoldChangeRhyth <- function(jdf){
  # Calculate log2 fold change between "rhythmic" and "non rhythmic" genes 
}

cossim <- function(x, y){
  return(x %*% y / sqrt(x%*%x * y%*%y))
}

LoopCor <- function(m, show.which=FALSE, input.vec1=NA){
  imin <- NA
  jmin <- NA
  jcor.min <- 2  # init because pearson cor is between -1 and 1. Use a number outside of this range.
  for (i in 1:ncol(m)){
    if (!is.na(input.vec1)){
      i <- input.vec1
    }
    vec1 <- m[, i]
    for (j in i:ncol(m)){
      if (j == i){
        next  # dont need to compare between same vec
      }
      vec2 <- m[, j]
      jcor <- cossim(vec1, vec2)
      if (jcor < jcor.min){
        jcor.min <- jcor
        imin <- i
        jmin <- j
      }
    }
    if (!is.na(input.vec1)){
      break
    }
  }
  if (show.which){
    jtissues <- c(colnames(m)[imin], colnames(m)[jmin])
    print(jtissues)
    return(list(tissues = jtissues,
                min.cor = jcor.min))
  } 
  return(jcor.min)
}

GetMinCor <- function(df){
  m <- acast(data = df, transcript ~ tissue, value.var = "norm_reads")
  jcor.min <- LoopCor(m)
  return(data.frame(min.cor = jcor.min))
}

Normalize <- function(x, pseudocount = 1){
  if (pseudocount > 0){
    x <- x + pseudocount
  }
  x.norm <- x / sum(x)
  return(x.norm)
}

AvgAcrossTranscripts <- function(df){
  ddply(df, .(transcript, gene, tissue), summarise, mean_reads = mean(reads))
}

NormalizeReads <- function(df){
  # Normalize reads across transcripts
  ddply(df, .(gene), transform, norm_reads = Normalize(mean_reads))
}

GetGeneNameFromAnnot <- function(annot){
  # from gene_name=RP23-271O17.1;transcript_id=ENSMUST00000193812 retrieve RP23-27017.1
  gene.str <- strsplit(as.character(annot), ';')[[1]][[1]]
  gene.str <- strsplit(gene.str, '=')[[1]][[2]]
  return(gene.str)
}

GetTranscriptIDFromAnnot <- function(annot){
  # from gene_name=RP23-271O17.1;transcript_id=ENSMUST00000193812 retrieve ENSMUST...
  transcript.str <- strsplit(as.character(annot), ';')[[1]][[2]]
  transcript.str <- strsplit(transcript.str, '=')[[1]][[2]]
  return(transcript.str)
}

rowMax <- function(df){
  # Return vector of maximums from df
  return(apply(df, 1, max))
}

GetTissuesAFE <- function(x){
  # Get tissues from column names from Adr_CT22 format
  substr(x, "_")[[1]][[1]]
}

TissueMapping <- function(cov.to.rnaseq = TRUE){
  # Tissue names from coverage are slightly different from
  # tissue names from RNASeq data. Create the mapping between 
  # coverage to rnaseq tissue names
  list("Adr" = "Adr",
       "Aor" = "Aorta",
       "BFat" = "BFAT",
       "Bstm" = "BS",
       "Cer" = "Cere",
       "Hrt" = "Heart",
       "Hyp" = "Hypo",
       "Kid" = "Kidney",
       "Liv" = "Liver",
       "Lun" = "Lung",
       "Mus" = "Mus",
       "WFat" = "WFAT")
}

ConvertTissuesToMatchRnaSeq <- function(tissues){
  # make tissue names look like RNASeq column names using TissueMapping
  tissue.map <- TissueMapping()
  sapply(tissues, function(x){
    tissue.map[[x]]
  })
}

GetExprsAsVector <- function(dat, genes, tissuetime){
  # Given a vector of rownames and column names, extract
  # its corresponding element in the dat. Return as
  # a vector.
  if (length(genes) != length(tissuetime)) print("Genes and tissuetime is not same length")
  lookups <- mapply(function(x, y){
    return(dat[x, y])
  }, x = genes, y = tissuetime)
  return(lookups)
}

GetLocationFromAnnotation <- function(bed, gene_name, transcript_id){
  # Given bed, gene_name and transcript_id, return the chromo, start, end
  annot <- paste0("gene_name=", gene_name, ";transcript_id=", transcript_id)
  sub <- subset(bed, annotations == annot)
  # return as UCSC-style
  return(paste0(sub$chromosome, ":", sub$start, "-", sub$end))
}

SubsetBed <- function(bed, gene_name, transcript){
  # Subset bed based on grepping annotations from gene name
  if (missing(transcript)){
    return(bed[grepl(gene_name, bed$annotations), ])  
  } else if (missing(gene_name)){  
    return(bed[grepl(transcript, bed$annotations), ])  
  }
}

NormalizeBySum <- function(x){
  # Normalize a vector by its sum
  return(x / sum(x))
}

GetFirst <- function(x){
  return(x[1])
}

ShannonEntropy <- function(x.vec, normalize=FALSE){
  if (normalize){
    # should sum to 1
    x.vec <- x.vec / sum(x.vec)
  }
  entropy <- 0
  for (x in x.vec){
    entropy <- entropy + x * log2(1 / x)
  }
  entropy <- entropy / log2(length(x.vec))
  return(entropy)
}

MulitpleStarts <- function(df, min_dist = 500){
  # Check if df has multiple exons, ddply from bed
  if (nrow(df) <= 1){
    return(data.frame(MultiStart = FALSE))
  }
  dist <- diff(c(min(df$start), max(df$start)))
  if (dist >= min_dist){
    return(data.frame(MultiStart = TRUE))
  } else{
    return(data.frame(MultiStart = FALSE))
  }
}

SubsetMinPval <- function(jdf){
  # take row with minimum pval
  jdf <- jdf[order(jdf$pval), ]
  return(jdf[1, ])
}

FitRhythNonRhyth <- function(jdf, log2transf = FALSE){
  # check it contains both Rhythmic and NotRhythmic elements
  if (length(unique(jdf$rhythmic.or.not)) != 2){
    return(data.frame(int = NA, coef = NA, pval = NA))
  }
  if (log2transf){
    jform <- log2(norm_reads) ~ rhythmic.or.not
  } else {
    jform <- norm_reads ~ rhythmic.or.not
  }
  jfit <- lm(data = jdf, formula = jform)
  # int <- summary(jfit)$mat["(Intercept)", "Coef"]
  # jcoef <- summary(jfit)$mat["rhythmic.or.notRhythmic", "Coef"]
  int <- coef(jfit)[[1]]
  jcoef <- coef(jfit)[[2]]
  pval <- summary(jfit)$coefficients["rhythmic.or.notRhythmic", "Pr(>|t|)"]
  return(data.frame(int = int, coef = jcoef, pval = pval))
}

FitPromoterUsageToAmplitude <- function(dat){
  #   jweights <- 1 / sqrt(dat$tpm_norm.var)
  #   fit <- lm(formula = tpm_norm.avg ~ relamp, data = dat, weights = jweights)
  fit <- lm(formula = tpm_norm.avg ~ relamp, data = dat)
  int <- coef(fit)[["(Intercept)"]]
  relamp <- coef(fit)[["relamp"]]
  tpm_norm.range <- diff(range(dat$tpm_norm.avg))
  relamp.range <- diff(range(dat$relamp))
  # handle if relamp is "NA"
  if (is.na(fit$coefficients[["relamp"]])){
    return(data.frame(int = NA, relamp = NA, pval = NA, tpm_norm.range = NA, relamp.range = NA))
  }
  pval <- summary(fit)$coefficient["relamp", "Pr(>|t|)"]
  return(data.frame(int = int, relamp = relamp, pval = pval, tpm_norm.range = tpm_norm.range, relamp.range = relamp.range))
}
