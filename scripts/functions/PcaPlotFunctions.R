# Jake Yeung
# PCA plotting functions
# Periodogram plotting functions
# Nov 5 2014
# 
PlotLoadings <- function(y, title="Plot title") {
  # Given vector from PCA, plot vector and color by tissue.
  # Maybe give fancy legend
  plot(y, main=title, col=rep(1:12, each=24), type='o') 
}
