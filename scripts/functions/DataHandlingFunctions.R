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
  print(dat[1:N, 1:N])
}