RunFisherDHS <- function(dat, S, jmodel.col = "peak.type", jshow.table = FALSE, jcutoff = 0.5){
  sitecounts.hash <- hash(as.character(dat$peak), dat$sitecount)
  S$sitecount <- sapply(S$peak, AssignSitecount, sitecounts.hash)
  test <- FisherTestSitecounts(dat = S, cutoff = jcutoff, sitecount.col = "sitecount", 
                               model.col = jmodel.col, show.table=jshow.table)
  return(test)
}

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
  if (nrow(N.table) != 2 | ncol(N.table) != 2){
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

RunFisherOnPromoters <- function(N.promoter, foreground.models, background.models, cutoffs, show.plot = TRUE){
  jtiss <- foreground.models
  N.sub <- subset(N.promoter, model %in% c(foreground.models, background.models))
  if (length(background.models) == 1){
    N.sub$model <- sapply(N.sub$model, function(m){
      if (!m %in% jtiss){
        return("Flat")
      } else {
        return("Rhyth")
      }
    })
  } else if (length(background.models) > 1){
    N.sub$model <- sapply(N.sub$model, function(m){
      if (m != jtiss){
        return("Flat")
      } else {
        return("Rhyth")
      }
    })
  } else {
    warning("Length of models must be natural number")
  }
  N.ftest.all <- data.frame()
  for (cutoff in cutoffs){
    N.ftest <- N.sub %>%
      group_by(motif) %>%
      do(FisherTestSitecounts(., cutoff))
    N.ftest$cutoff <- cutoff
    N.ftest.all <- rbind(N.ftest.all, as.data.frame(N.ftest))
  }
  N.ftest.sum <- N.ftest.all %>%
    group_by(motif) %>%
    summarise(odds.ratio = mean(odds.ratio), p.value = mean(p.value))
  N.ftest.sum <- N.ftest.all %>%
    group_by(motif) %>%
    summarise(odds.ratio = mean(odds.ratio), p.value = mean(p.value))
  if(show.plot){
    print(ggplot(N.ftest.sum, aes(y = -log10(p.value), x = odds.ratio, label = motif)) + geom_point() + geom_text() + theme_bw(24) + 
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                  panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position="bottom"))
  }
  return(N.ftest.sum)
}