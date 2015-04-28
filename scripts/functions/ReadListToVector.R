# ReadListToVector.R

ReadListToVector <- function(fname, skip_first_row=FALSE){
  # read file name, one column with no colnames as a vector.
  vec <- read.table(fname)
  vec <- as.character(unlist(vec))
  if (skip_first_row){
    return(vec[2:length(vec)])
  } else{
    return(vec)
  }
}