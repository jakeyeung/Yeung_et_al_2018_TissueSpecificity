OrderPhaseMatrix <- function(complex.mat, order.by = "Liver"){
  # Given complex.mat, order rows such that the average angles are
  # increasing.
  # 
  # Args:
  # complex.mat -> complex matrix
  # 
  # Returns:
  # complex.mat.ordered <- complex matrix, rows ordered by 
  # average angle across columns
#   avg.phases <- apply(complex.mat, 1, function(x){
#     # re.avg <- sum(cos(Re(x))) / length(x)
#     # im.avg <- sum(sin(Re(x))) / length(x)
#     # vec.avg <- complex(real = re.avg, imaginary = im.avg)
#     return(Arg(vec.avg))
#   })
#   print(avg.phases)
  tissue.phases <- Arg(complex.mat[, tissue])
  return(complex.mat[order(tissue.phases, decreasing = FALSE), , drop=FALSE])
}
