FisherTestSitecounts <- function(dat, cutoff, show.table=FALSE){
  dat$has.motif <- sapply(dat$sitecount, function(s){
    if (s > cutoff){
      return(TRUE)
    } else {
      return(FALSE)
    }
  })
  N.table <- table(dat$has.motif, unlist(dat$model))
  if (nrow(N.table) != ncol(N.table)){
    return(data.frame(odds.ratio = NA, p.value = NA))
  }
  test <- fisher.test(N.table)
  if (show.table){
    print(N.table)
    print(test)
  }
  return(data.frame(odds.ratio = test$estimate, p.value = test$p.value))
}