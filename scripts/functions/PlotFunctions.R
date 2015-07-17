library(ggplot2)
library(grid)

PlotAmpPhase <- function(dat){
  # Expect amp and phase in data frame column names.
  # label as gene names
  ggplot(data = dat, aes(x = amp, y = phase, label = gene)) + 
    geom_point() + 
    geom_text(aes(x = amp, y = phase, size = amp)) + 
    coord_polar(theta = "y") +
    xlab("Phase") +
    ylab("Amp") +
    scale_y_continuous(limits = c(0, 24), breaks = seq(2, 24, 2))
}

PlotAmpPhaseAllTissues <- function(dat, jtissue){
  # Expect amp and phase in data frame column names.
  # label as gene names
  if (!missing(jtissue)){
    dat <- subset(dat, tissue == jtissue)
  }
    ggplot(data = dat, aes(x = amp, y = phase, label = gene)) + 
      geom_point(size=0.5) + 
      geom_text(aes(x = amp, y = phase, size = amp)) + 
      coord_polar(theta = "y") +
      xlab("Phase") +
      ylab("Amp") +
      facet_wrap(~tissue) + 
      scale_y_continuous(limits = c(0, 24), breaks = seq(2, 24, 2)) 
}

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

PlotLoadings <- function(Loadings, title="Plot title", plot.colors, cex = 1) {
  # Given vector from PCA, plot vector and color by tissue.
  # Maybe give fancy legend
  if (missing(plot.colors)){
    plot.colors <- rep(1:12, each=24)
  }
  plot(Loadings, main=title, col=plot.colors, type='o',
       cex.axis = cex,
       cex.main = cex,
       cex.lab = cex)
}

PCbiplot <- function(PC, jtitle="Plot title", x="PC1", y="PC2") {
  # http://stackoverflow.com/questions/6578355/plotting-pca-biplot-with-ggplot2
  # PC being a prcomp object
  data <- data.frame(obsnames=row.names(PC$x), PC$x)
  data$labsize <- data[[x]] ^ 2 + data[[y]]^2
  plot <- ggplot(data, aes_string(x=x, y=y)) + geom_text(alpha=.3, aes(label=obsnames, size = labsize))
  plot <- plot + geom_hline(aes(0), size=.2) + geom_vline(aes(0), size=.2)
  datapc <- data.frame(varnames=rownames(PC$rotation), PC$rotation)
  mult <- min(
    (max(data[,y]) - min(data[,y])/(max(datapc[,y])-min(datapc[,y]))),
    (max(data[,x]) - min(data[,x])/(max(datapc[,x])-min(datapc[,x])))
  )
  datapc <- transform(datapc,
                      v1 = .7 * mult * (get(x)),
                      v2 = .7 * mult * (get(y))
  )
  plot <- plot + coord_equal() + geom_text(data=datapc, aes(x=v1, y=v2, label=varnames), size = 5, vjust=1, color="red")
  plot <- plot + geom_segment(data=datapc, aes(x=0, y=0, xend=v1, yend=v2), arrow=arrow(length=unit(0.2,"cm")), color="red")
  plot <- plot + ggtitle(jtitle)
  return(plot)
}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  # http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
  # Multiple plot function
  #
  # ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
  # - cols:   Number of columns in layout
  # - layout: A matrix specifying the layout. If present, 'cols' is ignored.
  #
  # If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
  # then plot 1 will go in the upper left, 2 will go in the upper right, and
  # 3 will go all the way across the bottom.
  #
  
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


# Add plot activities functions here --------------------------------------

PlotMeanActivitiesWithSE.singleintercept <- function(df){
  jgene <- unique(df$gene)
  ggplot(df,
         aes(x = tissue, y = exprs)) + 
    geom_line() + 
    geom_errorbar(aes(ymax = exprs + se, ymin = exprs - se)) +
    ylab("Activity") +
    ggtitle(jgene)
}

ConvertArgToPhase <- function(phase.rads, omega){
  # convert phase in radians to phase in time, using omega.
  # expect phase to be between -pi to pi, convert that
  # to 0 to 24 hours.
  
  # convert range from -pi to pi to 0 to 2pi
  phase.rads[which(phase.rads < 0)] <- phase.rads[which(phase.rads < 0)] + 2 * pi
  phase <- phase.rads / omega
  return(phase)
}

PlotComplex2 <- function(vec.complex, labels, omega = 2 * pi / 24, title = "My title", ylab = "Amplitude of activity", xlab = "Phase of activity (CT)"){
  # Convert complex to amplitude (2 * fourier amplitude) and phase, given omega.
  # then plot in polar coordinates
  # fourier amplitudes are half-amplitudes of the sine-wave
  # http://www.prosig.com/signal-processing/FourierAnalysis.pdf
  df <- data.frame(amp = Mod(vec.complex) * 2,
                   phase = ConvertArgToPhase(Arg(vec.complex), omega = omega),
                   label = labels)
  m <- ggplot(data = df, aes(x = amp, y = phase, label = label)) + 
    geom_point(size = 0.5) +
    coord_polar(theta = "y") + 
    xlab(xlab) +
    ylab(ylab) +
    geom_text(aes(x = amp, y = phase, size = amp), vjust = 0) +
    ggtitle(title) +
    scale_y_continuous(limits = c(0, 24), breaks = seq(2, 24, 2))
}

PlotEigengene <- function(svd.obj, comp, omega = 2 * pi / 24, rotate=TRUE){
  eigengene <- svd.obj$v[, comp]
  if (rotate){
    # rotate to phase of largest magnitude in sample of eigengene
    phase.reference <- Arg(eigengene[which(Mod(eigengene) == max(Mod(eigengene)))])
    rotate.factor <- complex(modulus = 1, argument = phase.reference)
    # rotate eigengene by -phase ref
    eigengene <- eigengene * Conj(rotate.factor)
  }
  v.plot <- PlotComplex2(eigengene, labels = rownames(svd.obj$v), omega = omega, title = paste("Right singular value", comp))
  print(v.plot)
}

PlotEigensamp <- function(svd.obj, comp, omega = 2 * pi / 24, rotate=TRUE){
  if (rotate){
    eigengene <- svd.obj$v[, comp]
    # rotate to phase of largest magnitude in sample of eigengene
    phase.reference <- Arg(eigengene[which(Mod(eigengene) == max(Mod(eigengene)))])
    rotate.factor <- complex(modulus = 1, argument = phase.reference)
  }
  eigensamp <- svd.obj$u[, comp]
  eigensamp <- eigensamp * rotate.factor
  u.plot <- PlotComplex2(eigensamp, labels = rownames(svd.obj$u), omega = omega, title = paste("Left singular value", comp))   
  print(u.plot)
}


# Heatmaps ----------------------------------------------------------------

PlotRelampHeatmap <- function(M, jtitle = "Plot Title"){
  library(gplots)
  my.palette <- colorRampPalette(c("black", "yellow"))(n = 300)
  # # (optional) defines the color breaks manually for a "skewed" color transition
  col.breaks = c(seq(0, 0.15, length=150),  # black
                 seq(0.151, 1, length=151))  # yellow
  heatmap.2(as.matrix(M), 
            col=my.palette, 
            breaks = col.breaks, 
            scale="none", 
            key=T, 
            keysize=1.5, 
            density.info = "density", 
            trace="none", 
            cexCol=0.9, 
            labRow=NA, 
            dendrogram = "both",
            main = jtitle)  
}
