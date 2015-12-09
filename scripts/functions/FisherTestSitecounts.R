FisherTestSitecounts <- function(dat, cutoff, sitecount.col, model.col, show.table=FALSE){
  # expects dat to have has.motif, sitecount
  if (missing(sitecount.col)){
    sitecount.col <- "sitecount"
  }
  if (missing(model.col)){
    model.col <- "model"
  }
  dat$has.motif <- sapply(dat[[sitecount.col]], function(s){
    if (s > cutoff){
      return(TRUE)
    } else {
      return(FALSE)
    }
  })
  N.table <- table(dat$has.motif, unlist(dat[[model.col]]))
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

SubsetAndFishers <- function(dat, jmodel, cutoffs){
  if (missing(cutoffs)){
    cutoffs <- seq(from = 0.4, to = 0.8, by = 0.1)
  }
  
  dat$model <- sapply(dat$model, function(m){
    if (!m %in% jmodel){
      return("Flat")
    } else {
      return("Rhyth")
    }
  })
  
  N.ftest.all <- data.frame()
  for (cutoff in cutoffs){
    print(cutoff)
    N.ftest <- dat %>%
      group_by(motif) %>%
      do(FisherTestSitecounts(., cutoff))
    N.ftest$cutoff <- cutoff
    N.ftest.all <- rbind(N.ftest.all, N.ftest)
  }
  
  N.ftest.sum <- N.ftest.all %>%
    group_by(motif) %>%
    summarise(odds.ratio = mean(odds.ratio), p.value = mean(p.value))
  return(N.ftest.sum)
}