# ReadListToVector.R

ReadListToVector <- function(fname){
  # read file name, one column with no colnames as a vector.
  vec <- read.table(fname)
  vec <- as.character(unlist(vec))
  return(vec)
}