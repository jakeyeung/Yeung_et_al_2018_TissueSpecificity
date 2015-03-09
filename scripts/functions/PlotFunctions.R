PlotComplex <- function(complex.matrix, gene.list, labels,
                        axis.min, axis.max, col="HSV",
                        main='Plot title', 
                        rotate=0,
                        add.text.plot=TRUE,
                        jpch=20,
                        threshold=0,
                        verbose=FALSE,
                        jcex=1){
  # Plot genes on complex plane
  # 
  # ARGS:
  # complex matrix Gene expression. Genes are rows, samples are columns.
  #   expect rownames in complex matrix whcih are gene names (and match genelist)
  # gene.list Optionally filter by gene list
  # colors Colors, if "HSV" compute HSV from angles using PhaseToHsv
  # axis.min: x and y axis min
  # axis.max: x and y axis max
  # main: title
  # rotate: rotate counter clockwise by how many radians? This gives start.angle
  # add.text.plot: add an optional text plot: useful if your gene list is huge.
  # jpch: pch plot symbols '.' for dot, 1 for circles, 20 for filled circles.
  # threshold: if many datapoints, show text label only if magnitude is greater than threshold
  # jcex: adjusting size of main, label and axis of plot. Useful for presentations.
  # 
  # requires wordcloud package for textplot function
  
  library(PhaseHSV)  # for colors around the circle
  
  if (add.text.plot){
    library(wordcloud)  # install.packages("wordcloud") 
  }
  
  # check if data is a matrix. If data frame, coerce to matrix.
  # this prevents errors: "Error non-numeric argument to function
  if (is.data.frame(complex.matrix)){
    complex.matrix <- as.matrix(complex.matrix)
  }
  
  if (missing(gene.list)){
    dat <- complex.matrix  
  } else {
    dat <- complex.matrix[gene.list, ]
  }
  if (missing(labels)){
    text.labels <- rownames(dat)  
  } else {
    text.labels <- labels
  }
  if (missing(axis.min) & missing(axis.max)){
    jmax <- max(Mod(complex.matrix))
    axis.min = -jmax
    axis.max = jmax
  }
  
  # rotate 
  rotation <- complex(modulus = 1, argument = rotate)
  dat <- dat * rotation
  
  plot.colors <- hsv(h=PhaseToHsv(Arg(dat), -pi, pi), s=1, v=1)
  
  plot(dat, col=plot.colors, 
       xlim=c(axis.min, axis.max), 
       ylim=c(axis.min, axis.max), 
       pch=jpch,
       main=main,
       xlab="Real",
       ylab="Complex",
       cex.main=jcex,
       cex.axis=jcex,
       cex.lab=jcex)
  abline(v=0)
  abline(h=0)
  
  # too many data points, only show labels for "large" datapoints
  filter.i <- which(Mod(dat) > threshold)
  text.labels[filter.i]
  text(dat[filter.i],
       labels=text.labels[filter.i], 
       pos=3)
  
  if (verbose){
    cat(paste0(text.labels[filter.i], collapse='", "'))
    # cat(paste0(text.labels[filter.i], collapse="\n"))
    cat("\n")
  }
  
  if (add.text.plot){
    if (is.null(dim(dat)) || nrow(dat) == 1) {
      warning("Data is not a matrix with >1 row. Not adding text plot")
    } else {
      textplot(Re(dat), Im(dat), text.labels, main=main,
               xlim=c(axis.min, axis.max),
               ylim=c(axis.min, axis.max),
               xlab="Real",
               ylab="Complex",
               cex=0.6,
               cex.main=jcex,
               cex.axis=jcex,
               cex.lab=jcex)
      abline(v=0)
      abline(h=0)  
    }
  }
}



PlotRnaMicroarrayFit <- function(tissue, gene, coeff.mat, array.exprs, rna.seq.exprs, rna.tissue){
  # Plotting function when fitting array with expression RNA Seq with microarray.
  # 
  # Args:
  # tissue: string representing tissue name in ARRAY e.g. "Liver"
  # gene: gene to plot
  # coeff.mat: contains intercept and slope for each tissue for every gene. Generated from
  # fit_array_with_exprs.rna.seq.vs.array.R
  # 
  # array.exprs: expression of array
  # 
  # rna.seq.exprs: expression in rnaseq
  # 
  # rna.tissue: string representing tissue name in rna.seq .eg. "Liv"
  # 
  # 
  
  # plot for one tissue only
  t.grep <- tissue
  t.grep.rnaseq <- rna.tissue
  
  intercept <- coeff.mat[gene, paste0(tissue, '_intercept')]
  slope <- coeff.mat[gene, paste0(tissue, '_slope')]
  
  x <- as.matrix(array.exprs[gene, which(grepl(t.grep, colnames(array.exprs)))])
  y <- as.matrix(rna.seq.exprs[gene, which(grepl(t.grep.rnaseq, colnames(rna.seq.exprs)))])
  plot(x, y, xlab='Array exprs',
       ylab='RNA exprs',
       main=paste(gene, 'Intercept=', signif(intercept, digits=3), 'Slope=', signif(slope, digits=3)))
  abline(intercept, slope, lty='dotted')
}

PlotComplexCircle <- function(complex.vector,
                              jlabels,
                              filter = 0,
                              period = 24,
                              rotate = 0,  # by hours
                              xlabel = "Magnitude",
                              ylabel = "Phase (hr)",
                              size = 24,
                              textsize = 6){
  library(ggplot2)
  
  theme_set(theme_grey(base_size = size))
  
  if (missing(jlabels)){
    no.label <- TRUE
  } else {
    no.label <- FALSE
  }
  
  jnames <- names(complex.vector)
  if (is.null(jnames)){
    jnames <- rep(NA, length(complex.vector))
  }
  complex.vector <- as.matrix(complex.vector)
  magnitude <- Mod(complex.vector)
  phase <- Arg(complex.vector)
  phase <- phase * (period / (2 * pi)) + rotate
  phase <- phase %% 24
  
  if (missing(jlabels)){
    jlabels <- jnames
  } else if (filter != 0){
    jlabels <- ifelse(magnitude < filter, "", jnames)
  } 
  
  dat <- data.frame(Magnitude = as.numeric(magnitude), Phase = as.numeric(phase), Labels = jlabels)
  jbreaks <- seq(length.out = period / 2, by = 2, to = period)
  if (!no.label){
    ggplot(data = dat, aes(x = Magnitude, y = Phase, label = Labels)) + 
      geom_point() + 
      geom_text(size = textsize) + 
      coord_polar(theta = "y") +
      xlab(xlabel) +
      ylab(ylabel) +
      scale_y_continuous(limits = c(0, 24), breaks = seq(2, 24, 2))
#   } else if (!no.label & filter > 0) {
#     ggplot(data = dat, aes(x = Magnitude, y = Phase, label = Labels)) + 
#       geom_point() + 
#       geom_text(size = textsize) + 
#       coord_polar(theta = "y") +
#       xlab(xlabel) +
#       ylab(ylabel) +
#       scale_y_continuous(limits = c(0, 24), breaks = seq(2, 24, 2))
  }  else if (no.label) {
    ggplot(data = dat, aes(x = Magnitude, y = Phase)) + 
      geom_point() + 
      coord_polar(theta = "y") +
      xlab(xlabel) +
      ylab(ylabel) +
      scale_y_continuous(limits = c(0, 24), breaks = seq(2, 24, 2))
  }
}

PlotFitDiagnostics <- function(array.exprs.full,
                               array.exprs.subset,
                               rna.seq.exprs, 
                               clockgenes, 
                               coeff.mat,
                               outpath, 
                               tissue, 
                               tissue.rna.seq,
                               allow.negs=FALSE){
  if (missing(tissue.rna.seq) == TRUE){
    tissue.rna.seq <- tissue
  }
    
  pdf(outpath)
  par(mfrow=c(2, 2))
  for (gene in clockgenes){
    slope <- coeff.mat[gene, paste0(tissue, "_slope")]
    intercept <- coeff.mat[gene, paste0(tissue, "_intercept")]
    
    # get microarray and rnaseq
    x <- array.exprs.full[gene, grepl(tissue, colnames(array.exprs))]
    y.rna.seq <- rna.seq.exprs[gene, grepl(tissue.rna.seq, colnames(rna.seq.exprs))]
    # get indices of samples with both microarray and rnaseq
    x.rna.seq.i <- colnames(array.exprs.subset[gene, grepl(tissue, colnames(array.exprs.subset))])  # subset only ones that have mRNA matched
    # convert to indices
    x.rna.seq.i <- names(x) %in% x.rna.seq.i  # logical True/False
    # get subset x, containing microarray and rnaseq
    x.subset <- x[x.rna.seq.i]
    
    plot.symbols <- sapply(x.rna.seq.i, function(true.false){
      if (true.false == TRUE){
        symbol <- 1
      } else {
        symbol <- 8
      }
      })
    y <- slope * x + intercept
    if (allow.negs == FALSE){
      y[which(y < 1)] <- 1
    }
    plot(2^x, 2^y, main=paste("o=samp w/ RNASeq+array", gene), pch=c(plot.symbols), xlab="Observed microarray (normal)", ylab="Predicted expression (normal)")
    abline(h=0)
    
    plot(x, y, main=paste("log2", gene, tissue), pch=c(plot.symbols), xlab="Observed microarray (log2)", ylab="Predicted expression (log2)")
    abline(h=0)
    
    # Plot raw on normal scale
    
    plot(2^x.subset, 2^y.rna.seq, main=paste("RNASeq vs Array: normal scale"), xlab="Microarray normal scale", ylab="DESeq normalized counts")
    
    # Plot raw on log2 scale
    plot(x.subset, y.rna.seq, 
         main=paste(gene, 'Intercept=', signif(intercept, digits=3), 
                    'Slope=', signif(slope, digits=3)), 
         xlab="Microarray log2 scale", 
         ylab="DESeq normalized counts (log2)") 
    abline(intercept, slope, lty='dotted')
    }
  dev.off() 
}

PlotBeforeAfter <- function(gene, array.before, array.after, rna.seq, y.max=14, convert.log2=FALSE, N.TISSUES=12){
  par(mfrow=c(3,1))
  if (convert.log2 == FALSE){
    array.before.gene <- as.numeric(array.before[gene, ])
    array.after.gene <- as.numeric(array.after[gene, ])
    rna.seq.gene <- as.numeric(rna.seq[gene, ])    
  } else if (convert.log2 == TRUE){
    array.before.gene <- log2(as.numeric(array.before[gene, ]))
    array.after.gene <- log2(as.numeric(array.after[gene, ]) + 1)
    rna.seq.gene <- log2(as.numeric(rna.seq[gene, ]) + 1)
  }

  # Array before adjustment
  plot(array.before.gene, main=paste(gene, 'log2 expression: array before adjustment'),
       col=rep(1:N.TISSUES, each=24), type='b', ylim=c(0, y.max), ylab="log2 exprs", 
       xlab=paste(tissue.names, collapse=" "))
  # Array after adjustment
  plot(array.after.gene, main=paste(gene, 'log2 exprs: array after adjustment'),
       col=rep(1:N.TISSUES, each=24), type='b', ylim=c(0, y.max), ylab="log2 exprs", 
       xlab=paste(tissue.names, collapse=" "))
  # RNA Seq
  plot(rna.seq.gene, main=paste(gene, 'log2 exprs: rnaseq'),
       col=rep(1:N.TISSUES, each=8), type='b', ylim=c(0, y.max), ylab="log2 exprs", 
       xlab=paste(tissue.names, collapse=" "))
  par(mfrow=c(1,1))
}

PlotAgainstRnaSeq <- function(gene, rna.seq, array.exprs.adjusted, 
                              common.samples, y.max=14){
  rna.seq.full <- matrix(NA, nrow=1, ncol=ncol(array.exprs.adjusted), 
                    dimnames=list(gene, colnames(array.exprs.adjusted)))
  rna.seq.full[gene, common.samples] <- as.matrix(rna.seq[gene, ])
  
  plot(seq(1:length(rna.seq.full)), rna.seq.full, main=paste(gene, 'black=rnaseq, red=array after adjust'),
       col=1, type='b', ylim=c(0, y.max), ylab="log2 exprs", 
       xlab=paste(tissue.names, collapse=" "))
  lines(as.matrix(array.exprs.adjusted[gene, ]), col=2, pch=22, type='o', cex=0.25)
}

GetFullR <- function(gene, rna.seq.exprs, common.samples){
  R.full <- matrix(0, nrow=1, ncol=ncol(array.exprs), 
                   dimnames=list(gene, colnames(array.exprs)))
  R.full[gene, common.samples] <- as.matrix(rna.seq.exprs[gene, common.samples])
  return(R.full)
}

GetUnobsObsSymbol <- function(all.samples, common.samples, unobs=8, obs=1){
  # create plot symbols. 8 = * = unobserved. 1 = o = observed.
  unobs.symbol <- unobs
  obs.symbol <- obs
  plot.symbols <- matrix(unobs.symbol, nrow=1, ncol=ncol(array.exprs),
                         dimnames=list(gene, colnames(array.exprs)))
  plot.symbols[gene, common.samples] <- obs.symbol
  return(plot.symbols)
}

PlotDiagnostics <- function(gene, array.exprs, rna.seq.exprs, 
                            common.samples, slope, int){
  # use full array.exprs, fill missing rna.seq.exprs with 0s and label with *
  # slope and int comes from coeff.mat
  # create R vs M full 288, R = 0 for "missing" values... for plotting
  
  # plot 2 by 1
  par(mfrow=c(2,1))
  
  R.full <- matrix(0, nrow=1, ncol=ncol(array.exprs), 
                   dimnames=list(gene, colnames(array.exprs)))
  R.full[gene, common.samples] <- as.matrix(rna.seq.exprs[gene, common.samples])
  M.full <- as.matrix(array.exprs[gene, ])
  # create M and R for fitting...
  R <- as.matrix(rna.seq.exprs[gene, common.samples])
  M <- as.matrix(array.exprs.subset.common.g[gene, ])
  # create plot symbols. 8 = * = unobserved. 1 = o = observed.
  unobs.symbol <- 8
  obs.symbol <- 1
  unobs.size <- 0.25
  obs.size <- 1
  plot.symbols <- matrix(unobs.symbol, nrow=1, ncol=ncol(array.exprs),
                         dimnames=list(gene, colnames(array.exprs)))
  plot.symbols[gene, common.samples] <- obs.symbol
  plot.cex <- sapply(plot.symbols, function(x){
    if(x == unobs.symbol){
      # 8 is unobserved 
      return(unobs.size)  # make half size
    } else {
      # 1 is observed
      return(obs.size)  # don't change size
    }
  })
  # int <- coeff.mat[gene, "intercept"]
  # slope <- coeff.mat[gene, "slope"]
  plot(M.full, R.full, main=paste(gene, "slope=", int, "int=", slope), 
       xlab="Microarray (log2)", 
       ylab="RNA-Seq (log2)",
       pch=plot.symbols,
       cex=plot.cex)
  abline(int, slope)
  # plot data on normal scale
  m.norm <- 2^M.full
  r.norm = 2^R.full - 1
  # draw its line in normal scale
  # convert log(y) = a * log(x) + b to normal scale
  f.r.norm <- function(slope, int, m) 2^(int) * m ^ slope - 1
  m.norm.predict <- seq(min(m.norm), max(m.norm), 10)
  r.norm.predict <- f.r.norm(slope, int, m.norm.predict)
  y.max <- max(c(r.norm, r.norm.predict))
  
  plot(m.norm, r.norm, main="norm. scale data. o=observed, *=unobserved",
       xlab="Microarray (normal scale)",
       ylab="RNA-seq (DESeq-normalized count",
       pch=plot.symbols,
       cex=plot.cex,
       ylim=c(0, y.max))
  lines(m.norm.predict, r.norm.predict)
  
  par(mfrow=c(1,1))
}

PlotArgsMatrix <- function(complex.mat, colors, main = "Title", jcex = 1){
  if (missing(colors)){
    hues <- seq(from=0, to=1, length.out=100)
    colors <- hsv(h=hues, s=1, v=1)
  }
  # --------- BEGIN: PLOT MATRIX OF PHASE ANGLES ----------------- # 
  # 
  # order by phase angle
  y <- 1:ncol(complex.mat)  # length of 12, genes are mix of these.
  x <- 1:nrow(complex.mat)  # length of top.N. samples are mix of these.
  par(mar=c(6.1, 5.1, 4.1, 2.1))
  image(x, y, Arg(complex.mat),
        col=colors,
        main=main, 
        axes=FALSE, xlab="", ylab="",
        cex.lab = jcex,
        cex.main = jcex,
        cex.axis = jcex)
  axis(1, at=x, labels=FALSE, tick=FALSE)
  axis(2, at=y, labels=FALSE, tick=FALSE)
  # Now draw the textual axis labels
  # gene labels
  if (jcex != 1){
    gene.cex <- 1
  } else {
    gene.cex <- 0.5
  }
  text(x, par("usr")[3] - 1.5,
       labels=rownames(complex.mat),
       srt=90, 
       pos=4,
       offset=0,
       xpd=TRUE, 
       cex=gene.cex) 
  # sample labels
  text(par("usr")[1] - .3, y, 
       labels = colnames(complex.mat), 
       srt=0, 
       pos=2, 
       offset=0,
       xpd=TRUE, 
       cex=1.5) 
  # ---------- END: PLOT MATRIX OF PHASE ANGLES ------------------- # 
}