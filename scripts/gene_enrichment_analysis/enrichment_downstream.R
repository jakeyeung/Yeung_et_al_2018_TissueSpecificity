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

dir.create("plots/GO_analysis")
pdf("plots/GO_analysis/tissue_modules.pdf")

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# liver WT KO -------------------------------------------------------------


jmod <- "Liver_SV129,Liver_BmalKO"
load(paste0("Robjs/GO_analysis/model", jmod, ".Robj"), v=T); enrichment.liverWTKO <- CopyZT0ToZT24(enrichment); rm(enrichment)

ignore.GOs <- c("GO:0042254")  # ribi
ignore.GOs <- c("GO:0043434")  # peptide hrmone
liverWTKOsub <- subset(enrichment.liverWTKO, !GO.ID %in% ignore.GOs)
amp.max <- max(liverWTKOsub$minuslogpval)
ggplot(liverWTKOsub, 
       aes(x = tmid, y = minuslogpval, fill = Term)) + 
  geom_polygon(alpha = 0.3) + 
  coord_polar(theta = "x") + scale_x_continuous(limits = c(0, 24), breaks = seq(6, 24, 6)) + 
  scale_fill_manual(values = cbPalette) + 
  theme_bw() +
  ggtitle(jmod) + 
  geom_hline(yintercept = seq(0, amp.max, length.out = 2), colour = "grey50", size = 0.2, linetype = "dashed") +
  geom_vline(xintercept = seq(6, 24, by = 6), colour = "grey50", size = 0.2, linetype = "solid") +
  theme(panel.grid.major = element_line(size = 0.5, colour = "grey"), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position="bottom",
        panel.border = element_blank(),
        legend.key = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank())

ggplot(liverWTKOsub, 
       aes(x = tmid, y = minuslogpval, colour = Term)) + 
  geom_line() + 
  coord_polar(theta = "x") + scale_x_continuous(limits = c(0, 24), breaks = seq(6, 24, 6)) + 
  scale_fill_manual(values = cbPalette) + 
  theme_bw() +
  ggtitle(jmod) + 
  geom_hline(yintercept = seq(0, amp.max, length.out = 2), colour = "grey50", size = 0.2, linetype = "dashed") +
  geom_vline(xintercept = seq(6, 24, by = 6), colour = "grey50", size = 0.2, linetype = "solid") +
  theme(panel.grid.major = element_line(size = 0.5, colour = "grey"), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position="bottom",
        panel.border = element_blank(),
        legend.key = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank()) + 
  facet_wrap(~Term)


# Liver WT ----------------------------------------------------------------

jmod <- "Liver_SV129"
load(paste0("Robjs/GO_analysis/model", jmod, ".Robj"), v=T); enrichment.liverWT <- CopyZT0ToZT24(enrichment); rm(enrichment)
liverWTsub <- subset(enrichment.liverWT, !GO.ID %in% c("GO:0006633", "GO:0006511", "GO:0055114", "GO:0070542", "GO:0007584", "GO:0009056"))
amp.max <- max(liverWTsub$minuslogpval)
ggplot(liverWTsub, 
       aes(x = tmid, y = minuslogpval, fill = Term)) + 
  geom_polygon(alpha = 0.3) + 
  coord_polar(theta = "x") + scale_x_continuous(limits = c(0, 24), breaks = seq(6, 24, 6)) + 
  scale_fill_manual(values = cbPalette) + 
  theme_bw() +
  ggtitle(jmod) + 
  geom_hline(yintercept = seq(0, amp.max, length.out = 2), colour = "grey50", size = 0.2, linetype = "dashed") +
  geom_vline(xintercept = seq(6, 24, by = 6), colour = "grey50", size = 0.2, linetype = "solid") +
  theme(panel.grid.major = element_line(size = 0.5, colour = "grey"), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position="bottom",
        panel.border = element_blank(),
        legend.key = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank())

ggplot(liverWTsub, 
       aes(x = tmid, y = minuslogpval, colour = Term)) + 
  geom_line() + 
  coord_polar(theta = "x") + scale_x_continuous(limits = c(0, 24), breaks = seq(6, 24, 6)) + 
  scale_fill_manual(values = cbPalette) + 
  theme_bw() +
  ggtitle(jmod) + 
  geom_hline(yintercept = seq(0, amp.max, length.out = 2), colour = "grey50", size = 0.2, linetype = "dashed") +
  geom_vline(xintercept = seq(6, 24, by = 6), colour = "grey50", size = 0.2, linetype = "solid") +
  theme(panel.grid.major = element_line(size = 0.5, colour = "grey"), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position="bottom",
        panel.border = element_blank(),
        legend.key = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank()) + 
  facet_wrap(~Term)


# Kidney WT ---------------------------------------------------------------

jmod <- "Kidney_SV129"
load(paste0("Robjs/GO_analysis/model", jmod, ".Robj"), v=T); enrichment.kidneyWT <- CopyZT0ToZT24(enrichment); rm(enrichment)
kidneysub <- subset(enrichment.kidneyWT, !GO.ID %in% c("GO:0006633", "GO:0006511", "GO:0055114", "GO:0070542", "GO:0007584", "GO:0009056"))
ampmax <- max(kidneysub$minuslogpval)
ggplot(kidneysub, 
       aes(x = tmid, y = minuslogpval, fill = Term)) + 
  geom_polygon(alpha = 0.3) + 
  coord_polar(theta = "x") + scale_x_continuous(limits = c(0, 24), breaks = seq(6, 24, 6)) +
  scale_fill_manual(values = cbPalette) +
  theme_bw()  + 
  theme(legend.position = "bottom") + 
  ggtitle(jmod)

ggplot(kidneysub, aes(x = tmid, y = minuslogpval, colour = Term)) + 
  geom_line() + 
  coord_polar(theta = "x") + scale_x_continuous(limits = c(0, 24), breaks = seq(6, 24, 6)) + 
  scale_colour_manual(values = cbPalette) + 
  theme_bw()  + 
  theme(legend.position = "bottom") + 
  facet_wrap(~Term) + 
  ggtitle(jmod)


# Liver WTKO and Liver WT -------------------------------------------------

jmod <- "Liver_SV129,Liver_BmalKO-Liver_SV129"
load(paste0("Robjs/GO_analysis/model", jmod, ".Robj"), v=T); enrichment <- CopyZT0ToZT24(enrichment)
kidneysub <- subset(enrichment, !GO.ID %in% c("GO:0006633", "GO:0006511", "GO:0055114", "GO:0070542", "GO:0007584", "GO:0009056"))
ampmax <- max(kidneysub$minuslogpval)
ggplot(kidneysub, 
       aes(x = tmid, y = minuslogpval, fill = Term)) + 
  geom_polygon(alpha = 0.3) + 
  coord_polar(theta = "x") + scale_x_continuous(limits = c(0, 24), breaks = seq(6, 24, 6)) +
  # scale_fill_manual(values = cbPalette) +
  theme_bw()  + 
  theme(legend.position = "bottom") + 
  ggtitle(jmod)

ggplot(kidneysub, aes(x = tmid, y = minuslogpval, colour = Term)) + 
  geom_line() + 
  coord_polar(theta = "x") + scale_x_continuous(limits = c(0, 24), breaks = seq(6, 24, 6)) + 
  # scale_colour_manual(values = cbPalette) + 
  theme_bw()  + 
  theme(legend.position = "bottom") + 
  facet_wrap(~Term) + 
  ggtitle(jmod)


dev.off()