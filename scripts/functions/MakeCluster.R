MakeCluster <- function(ncores, rscript_args){
  library(doParallel)
  # Make cluster pass no args or else it will take FOREVER for cl to be made.
  # Ref: http://stackoverflow.com/questions/23540129/what-affects-the-time-to-create-a-cluster-using-the-parallel-package
  if (missing(ncores)){
    ncores <- detectCores()
  }
  if (missing(rscript_args)){
    rscript_args = c("--no-init-file", "--no-site-file", "--no-environ")
  }
  cl <- makeCluster(ncores, rscript_args=rscript_args)
  registerDoParallel(cl)
}