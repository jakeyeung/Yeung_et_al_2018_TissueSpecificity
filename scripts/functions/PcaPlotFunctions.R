# Jake Yeung
# PCA plotting functions
# Periodogram plotting functions
# Nov 5 2014
# 
PlotLoadings <- function(Loadings, title="Plot title", cex = 1) {
  # Given vector from PCA, plot vector and color by tissue.
  # Maybe give fancy legend
  plot(Loadings, main=title, col=rep(1:12, each=24), type='o',
       cex.axis = cex,
       cex.main = cex,
       cex.lab = cex)
}
