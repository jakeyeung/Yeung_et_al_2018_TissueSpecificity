# Jake Yeung
# quality.control.normalization.R
# Dec 2 2014

# Define dirs -------------------------------------------------------------

# define dirs
data.dir <- "data"
normalized.array.fname <- file.path(data.dir, "array.adjusted.lm.log2.txt")


# Load file ---------------------------------------------------------------

normalized.array <- read.table(normalized.array.fname)


# How many have negative values? ------------------------------------------

negs <- apply(normalized.array, 1, function(x){
  if (min(x) < 0){
    return(1)
  } else {
    return(0)
  }
})

problem.genes <- names(negs[which(negs == 1)])

# pick random
set.seed(0)
lucky.genes <- sample(problem.genes, 30)


# All negatives -> correct to 0. ------------------------------------------


