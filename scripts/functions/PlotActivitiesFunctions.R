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