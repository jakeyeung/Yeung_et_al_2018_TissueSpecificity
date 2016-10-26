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

ignore.GOs <- c("GO:0042254")  # ribi
ignore.GOs <- c("GO:0043434")  # peptide hrmone
liverWTKOsub <- subset(enrichment.liverWTKO, !GO.ID %in% ignore.GOs)
ggplot(liverWTKOsub, 
       aes(x = tmid, y = minuslogpval, fill = Term)) + 
  geom_polygon(alpha = 0.3) + 
  coord_polar(theta = "x") + scale_x_continuous(limits = c(0, 24), breaks = seq(6, 24, 6)) + 
  scale_fill_manual(values = cbPalette) + 
  theme_bw() +
  theme(legend.position = "bottom")

ggplot(liverWTKOsub, 
       aes(x = tmid, y = minuslogpval, colour = Term)) + 
  geom_line() + 
  coord_polar(theta = "x") + scale_x_continuous(limits = c(0, 24), breaks = seq(6, 24, 6)) + 
  scale_fill_manual(values = cbPalette) + 
  theme_bw() +
  theme(legend.position = "bottom")

load("Robjs/GO_analysis/modelLiver_SV129.Robj", v=T); enrichment.liverWT <- CopyZT0ToZT24(enrichment); rm(enrichment)


liverWTsub <- subset(enrichment.liverWT, !GO.ID %in% c("GO:0006633", "GO:0006511", "GO:0055114", "GO:0070542", "GO:0007584", "GO:0009056"))
ggplot(liverWTsub, 
       aes(x = tmid, y = minuslogpval, fill = Term)) + 
  geom_polygon(alpha = 0.3) + 
  coord_polar(theta = "x") + scale_x_continuous(limits = c(0, 24), breaks = seq(6, 24, 6)) + 
  scale_fill_manual(values = cbPalette) + 
  theme_bw()  + 
  theme(legend.position = "bottom")

ggplot(liverWTsub, aes(x = tmid, y = minuslogpval, colour = Term)) + 
  geom_line(alpha = 0.3) + 
  coord_polar(theta = "x") + scale_x_continuous(limits = c(0, 24), breaks = seq(6, 24, 6)) + 
  scale_colour_manual(values = cbPalette) + 
  theme_bw()  + 
  theme(legend.position = "bottom")

# create a single color scheme  
