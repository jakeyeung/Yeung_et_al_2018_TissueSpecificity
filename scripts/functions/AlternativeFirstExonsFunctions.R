# Functions: kallisto -----------------------------------------------------

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
