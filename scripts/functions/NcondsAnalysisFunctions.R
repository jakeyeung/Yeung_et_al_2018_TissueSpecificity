GetTissuesFromCoefFit <- function(fit.coef.names, tissue.colname = "tissue"){
  # Get tissues involved by grepping everything after tissue.colname
  match.names.i <- grepl(tissue.colname, fit.coef.names) & !grepl(":", fit.coef.names)
  match.names <- fit.coef.names[match.names.i]
  tissues <- sapply(match.names, function(jname) strsplit(jname, split = tissue.colname)[[1]][[2]], USE.NAMES = FALSE)
  return(tissues)
}

ExtraParamsFromFit <- function(fit.coef, tissues){
  # From coefficients of fit, extract parameters
  # colnames should contain tissue, RNA-Seq intercept, phase, amplitude, and BIC
  cnames <- tissues
  # TODO
}

SubsetByMaxBicWeight <- function(dat){
  max.i <- which.max(dat$bicweight)
  return(dat[max.i, ])
}

PlotSvdFromGeneList <- function(dat.complex, gene.list, jcomp = 1, jlabel.n = 25, jvalue.var = "exprs.transformed", constant.amp = 4){
  # Plot SVD from gene list
  jlayout <- matrix(c(1, 2, 3, 4), 2, 2, byrow = TRUE)
  s.custom <- SvdOnComplex(subset(dat.complex, gene %in% gene.list), value.var = jvalue.var)
  eigens.custom <- GetEigens(s.custom, period = 24, comp = jcomp, label.n = jlabel.n, eigenval = TRUE, adj.mag = TRUE, constant.amp = constant.amp)
  eigens.custom2 <- GetEigens(s.custom, period = 24, comp = (jcomp + 1), label.n = jlabel.n, eigenval = TRUE, adj.mag = TRUE, constant.amp = constant.amp)
  multiplot(eigens.custom$u.plot, eigens.custom$v.plot, eigens.custom2$u.plot, eigens.custom2$v.plot, layout = jlayout)
  return(eigens.custom)
}

GetTopGenesFromSvd <- function(eigens.obj, top.n = 25){
  return(names(sort(eigens.obj$eigensamp, decreasing = TRUE)[1:top.n]))
}

PCbiplot <- function(PC, jsizes, x="PC1", y="PC2", max.size = 5) {
  if (missing(jsizes)) jsizes <- 1
  # PC being a prcomp object
  dat <- data.frame(obsnames=row.names(PC$x), size = jsizes, PC$x)
  plot <- ggplot(dat, aes_string(x=x, y=y)) + geom_text(alpha=.8, aes(label=obsnames, size = size)) + scale_size(range = c(0, max.size))
  plot <- plot + geom_hline(aes(0), size=.2) + geom_vline(aes(0), size=.2)
  datapc <- data.frame(varnames=rownames(PC$rotation), PC$rotation)
  mult <- min(
    (max(dat[,y]) - min(dat[,y])/(max(datapc[,y])-min(datapc[,y]))),
    (max(dat[,x]) - min(dat[,x])/(max(datapc[,x])-min(datapc[,x])))
  )
  datapc <- transform(datapc,
                      v1 = .7 * mult * (get(x)),
                      v2 = .7 * mult * (get(y))
  )
  plot <- plot + coord_equal() + geom_text(data=datapc, aes(x=v1, y=v2, label=varnames), size = 5, vjust=1, color="red")
  plot <- plot + geom_segment(data=datapc, aes(x=0, y=0, xend=v1, yend=v2), arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color="red")
  return(plot)
}

PlotTissuePcaFromGenelist <- function(dat.means, genelist, fits.comb){
  # dat.means: gene, tissue, exprs.mean long format
  # genelist: list of genes
  dat.mat.comb <- dcast(subset(dat.means, gene %in% genelist), formula = gene ~ tissue, value.var = "exprs.mean")
  rownames(dat.mat.comb) <- dat.mat.comb$gene; dat.mat.comb$gene <- NULL
  
  # center
  dat.mat.comb <- t(scale(t(dat.mat.comb), center = TRUE, scale = FALSE))
  dat.mat.pca <- prcomp(dat.mat.comb, center = FALSE, scale. = FALSE)
  
  #   biplot(dat.mat.pca)
  jsizes <- fits.comb$bicweight[order(fits.comb$gene)]
  print(PCbiplot(dat.mat.pca, jsizes))
}