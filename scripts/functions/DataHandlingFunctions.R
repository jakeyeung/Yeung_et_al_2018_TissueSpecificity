# Jake Yeung
# DataHandlingFunctions.R
# November 10 2014
# Simple functions for easier manipulating of data

Peek <- function(dat, N=5){
  # Peek at first N rows and N columns of data
  # 
  # ARGS
  #   dat: dataframe or matrix
  # 
  # OUTPUT
  #   print to user the first N rows and N columns of data
  if (N > nrow(dat)){
    N.row <- nrow(dat)
  } else {
    N.row <- N
  }
  if (N > ncol(dat)){
    N.col <- ncol(dat)
  } else {
    N.col <- N
  }
  print(dat[1:N.row, 1:N.col])
}