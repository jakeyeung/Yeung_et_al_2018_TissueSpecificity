ShannonEntropy <- function(p){
  s <- p * log2(1 / p)
  return(sum(s))
}