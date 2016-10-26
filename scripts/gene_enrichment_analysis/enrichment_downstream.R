# 2016-10-26
# after doing sliding window GO analysis: plot data

library(ggplot2)

CopyZT0ToZT24 <- function(enrichment, twindow = 6){
  # complete the circle
  enrichment$tmid <- enrichment$tstart + twindow / 2
  enrichment$tmid <- sapply(enrichment$tmid, function(tmid) ifelse(tmid >= 24, tmid - 24, tmid))
  
  enrichment.24 <- subset(enrichment, tmid == 0)
  enrichment.24$tmid <- 24
  enrichment <- bind_rows(enrichment, enrichment.24)
  return(enrichment)
}

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

load("Robjs/GO_analysis/modelLiver_SV129,Liver_BmalKO.Robj", v=T); enrichment.liverWTKO <- CopyZT0ToZT24(enrichment); rm(enrichment)
ggplot(enrichment.liverWTKO, aes(x = tmid, y = minuslogpval, colour = Term)) + geom_line() + coord_polar(theta = "x") + scale_x_continuous(limits = c(0, 24), breaks = seq(6, 24, 6)) + 
  scale_colour_manual(values = cbPalette)

load("Robjs/GO_analysis/modelLiver_SV129.Robj", v=T); enrichment.liverWT <- CopyZT0ToZT24(enrichment); rm(enrichment)
ggplot(subset(enrichment.liverWT, !GO.ID %in% c("GO:0006633", "GO:0006511")), aes(x = tmid, y = minuslogpval, colour = Term)) + 
         geom_line() + coord_polar(theta = "x") + scale_x_continuous(limits = c(0, 24), breaks = seq(6, 24, 6)) + 
  scale_colour_manual(values = cbPalette)
  
