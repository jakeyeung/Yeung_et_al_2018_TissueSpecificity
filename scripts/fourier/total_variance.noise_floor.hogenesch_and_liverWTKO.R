# 2017-02-09
# Jake Yeung
# assess noise floor in total variance
# for both Hogenesch and Liver WTKO

rm(list=ls())

library(ggplot2)
library(dplyr)

remove.wfat <- TRUE
barwidth <- 0.7
ymax <- 625
ymin <- 0
yinc <- 200
jylab <- "Variance [log2]^2"
jltype <- "dotted"

# uncomment out the ds you want
ds <- "Hogenesch"
ds <- "LiverWTKO"

# Functions ---------------------------------------------------------------

source("scripts/functions/FourierFunctions.R")

MakeTissueGeno <- function(dat){
  dat$tissue.tmp <- dat$tissue
  dat$tissue_full <- as.character(dat$tissue.tmp)
  dat$tissue <- sapply(dat$tissue_full, function(x) strsplit(x, "_")[[1]][[1]])
  dat$tissue <- factor(dat$tissue, levels = c("Liver", "Kidney"))
  dat$geno <- AddGeno(dat, dat$tissue_full, return.df = FALSE)
  dat$tissue.tmp <- NULL
  dat$tissue_full <- NULL
  return(dat)
}

AddGeno <- function(dat, tissue, return.df = FALSE){
  dat$geno <- sapply(as.character(tissue), function(x) strsplit(x, "_")[[1]][[2]])
  dat$geno <- RenameGeno(dat$geno, return.as.factor = TRUE)
  if (return.df){
    return(dat) 
  } else {
    return(dat$geno)
  }
}

RenameGeno <- function(x, return.as.factor = TRUE, ko.str = "Bmal1KO"){
  x <- gsub(pattern = "SV129", replacement = "WT", x)
  x <- gsub(pattern = "BmalKO", replacement = ko.str, x)
  if (return.as.factor){
    x <- factor(x, levels = c("WT", ko.str))
  }
  return(x)
}


UpdateStartEnd <- function(dat){
  # make x coordinates match for collapsed data
  x1s <- dat$start[[1]]
  x1e <- (dat$start[[1]] + dat$start[[1]]) / 2
  x2s <- x1e
  x2e <- dat$end[[1]]
  x3s <- dat$start[[2]]
  x3e <- (dat$start[[2]] + dat$end[[2]]) / 2
  x4s <- x3e
  x4e <- dat$end[[2]]
  dat$start <- c(x1s, x2s, x3s, x4s)
  dat$end <- c(x1e, x2e, x3e, x4e)
  return(dat)
}

# Select dataset ----------------------------------------------------------


if (ds == "Hogenesch"){
  load("/home/yeung/projects/tissue-specificity/Robjs/dat.complex.all_T.rbinded.Robj", v=T)
} else if (ds == "LiverWTKO"){
  load("/home/yeung/projects/tissue-specificity/Robjs/liver_kidney_atger_nestle/dat.complex.all_T.bugfixed.Robj", v=T)
} else {
  warning(paste("ds must be Hogenesch or LiverWTKO", "ds=", ds))
}
print(paste("Loaded data for ds=", ds))

pdf(paste0("/home/yeung/projects/tissue-specificity/plots/fourier_analysis_with_noise_floor/fourier_analysis_noise_floor.", ds, ".pdf"))

# Figure 1B Temporal variation using Fourier analysis ------------------

if (remove.wfat){
  dat.complex.all_T <- subset(dat.complex.all_T, tissue != "WFAT")
}

dat.var.s <- dat.complex.all_T %>%
  group_by(period, tissue) %>%
  summarise(sum_sqr_mod = sum(Mod(exprs.transformed) ^ 2) * 2) %>%  # * 2 to consider symmetrical frequencies
  mutate(period.factor = signif(period, digits = 3))
dat.var.s$period.factor <- factor(dat.var.s$period.factor, 
                                  levels = sort(unique(dat.var.s$period.factor), decreasing = TRUE))
dat.var.all <- dat.var.s %>%
  group_by(tissue) %>%
  summarise(sum_sqr_mod.total = sum(sum_sqr_mod) * 2)  # * 2 considers symmetrical frequencies
dat.var.all$tissue <- factor(dat.var.all$tissue,
                             levels = dat.var.all$tissue[order(dat.var.all$sum_sqr_mod.total, decreasing = TRUE)])

# period.factor: everything not 24 and 12 should be called Other
dat.var.s$period.factor.cond <- sapply(dat.var.s$period.factor, function(f){
  f <- as.character(f)
  if (f == "24" | f == "12"){
    return(as.factor(f))
  }
  else{
    return(as.factor("Other"))
  }
})
dat.var.s <- dat.var.s %>%
  group_by(tissue) %>%
  arrange(period.factor.cond)

dat.var.s1_adj <- dat.var.s %>%
  group_by(tissue) %>%
  mutate(s1_normalized = sum_sqr_mod / sum(sum_sqr_mod))

dat.var.s1_adj.24 <- subset(dat.var.s1_adj, period.factor.cond == "24")
dat.var.s1_adj.24$tissue <- factor(dat.var.s1_adj.24$tissue,
                                   levels = dat.var.s1_adj.24$tissue[order(dat.var.s1_adj.24$s1_normalized, decreasing = TRUE)])

dat.var.s1_adj.12 <- subset(dat.var.s1_adj, period.factor.cond == "12")
dat.var.s1_adj.12$tissue <- factor(dat.var.s1_adj.12$tissue,
                                   levels = dat.var.s1_adj.12$tissue[order(dat.var.s1_adj.12$s1_normalized, decreasing = TRUE)])

# for ordering the facet_wrap plot across tissues by total variance
dat.var.s$tissue <- factor(dat.var.s$tissue,
                           levels = dat.var.all$tissue[order(dat.var.all$sum_sqr_mod.total, decreasing = TRUE)])

# for ordering the facet_wrap plot across tissues by 24h variance
dat.var.s$tissue <- factor(dat.var.s$tissue,
                           levels = dat.var.s1_adj.24$tissue[order(dat.var.s1_adj.24$sum_sqr_mod, decreasing = TRUE)])
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")
barwidth=0.7
plot.spectral.power <- ggplot() + 
  geom_bar(data=subset(dat.var.s, period.factor.cond != "Other"), mapping=aes(x = tissue, y = sum_sqr_mod, fill = period.factor.cond),
           stat="identity", position=position_dodge(), colour="black", width = barwidth) +
  theme_bw(24) + theme(axis.text.x=element_text(angle=90,vjust = 0, hjust = 1), 
                       panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank(), 
                       panel.background = element_blank(), 
                       axis.line = element_line(colour = "black"),
                       aspect.ratio = 0.5, 
                       legend.position="bottom") + xlab("") + ylab(jylab) +
  scale_y_continuous(breaks = c(seq(ymin, ymax, by = yinc)), limits = c(ymin, ymax)) + 
  scale_fill_manual(name = "Period [h]", drop=TRUE, values=cbPalette)
# add noise floor
periods <- sort(unique(dat.var.s$period))
noise.components <- periods[which(24 %% periods != 0)]

# mean for all tissues
dat.var.s.noisefloor <- subset(dat.var.s, period %in% noise.components) %>%
  group_by(period.factor) %>%
  summarise(sum_sqr_mod.mean = mean(sum_sqr_mod), 
            sum_sqr_mod.var = var(sum_sqr_mod), 
            sum_sqr_mod.min = min(sum_sqr_mod), 
            sum_sqr_mod.max = max(sum_sqr_mod))
# tissue by tissue
dat.var.s.noisefloor.bytiss <- subset(dat.var.s, period %in% noise.components) %>%
  group_by(tissue) %>%
  summarise(sum_sqr_mod.mean = mean(sum_sqr_mod))
lstart <- 0.65
dat.var.s.noisefloor.bytiss$start <- seq(lstart, by = 1, length.out = nrow(dat.var.s.noisefloor.bytiss))
dat.var.s.noisefloor.bytiss$end <- seq(lstart + barwidth, by = 1, length.out = nrow(dat.var.s.noisefloor.bytiss))

# add noise floor that is tissuewide
noise.floor <- mean(dat.var.s.noisefloor$sum_sqr_mod.mean)
plot.spectral.power.noiseflr <- plot.spectral.power + geom_hline(aes(yintercept = noise.floor), linetype=jltype)

# add tissue-specific noise floor
plot.spectral.power.noiseflr.bytiss <- plot.spectral.power + geom_segment(data = dat.var.s.noisefloor.bytiss, mapping=aes(x=start, xend=end, y=sum_sqr_mod.mean, yend=sum_sqr_mod.mean), linetype = jltype)


# plot Fraction of Total Variance
dat.var.s1_adj$tissue <- factor(dat.var.s1_adj$tissue,
                                levels = dat.var.s1_adj.24$tissue[order(dat.var.s1_adj.24$s1_normalized, decreasing = TRUE)])
# add tissue-specific noise floor
dat.var.s1_adj.noisefloor.bytiss <- subset(dat.var.s1_adj, period %in% noise.components) %>%
  group_by(tissue) %>%
  summarise(s1_normalized = mean(s1_normalized))
dat.var.s1_adj.noisefloor.bytiss$start <- seq(lstart, by = 1, length.out = nrow(dat.var.s1_adj.noisefloor.bytiss))
dat.var.s1_adj.noisefloor.bytiss$end <- seq(lstart + barwidth, by = 1, length.out = nrow(dat.var.s1_adj.noisefloor.bytiss))

plot.normalized.spectral.power <- ggplot() + 
  geom_bar(data=subset(dat.var.s1_adj, period.factor.cond != "Other"), mapping=aes(x = tissue, y = s1_normalized, fill = period.factor.cond),
           stat="identity", position=position_dodge(), colour="black", width = barwidth) +
  theme_bw(24) + theme(axis.text.x=element_text(angle=90,vjust = 0, hjust = 1), 
                       panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank(), 
                       panel.background = element_blank(), 
                       axis.line = element_line(colour = "black"),
                       aspect.ratio = 0.5, 
                       legend.position="bottom") + xlab("") + ylab("Fraction of Total Variance") +
  scale_fill_manual(name = "Period [h]", drop=TRUE, values=cbPalette)
plot.normalized.spectral.power <- plot.normalized.spectral.power + geom_segment(data = dat.var.s1_adj.noisefloor.bytiss, mapping=aes(x=start, xend=end, y=s1_normalized, yend = s1_normalized), linetype = jltype)
### END FOURIER ANALYSIS ###

print(plot.spectral.power)
print(plot.spectral.power.noiseflr)
print(plot.spectral.power.noiseflr.bytiss)
print(plot.normalized.spectral.power)


# Compare WT and KO for Liver and Kidney ----------------------------------

if (ds == "LiverWTKO"){
  dat.var.s2 <- MakeTissueGeno(dat.var.s)
  dat.var.s.dotted <- AddGeno(dat.var.s, dat.var.s$tissue, return.df = TRUE)
  dat.var.s.dotted$tissue <- RenameGeno(dat.var.s.dotted$tissue, return.as.factor = FALSE, ko.str = "KO")  # ordering reversed
  dat.var.s.dotted$tissue <- gsub(pattern = "_", replacement = " ", dat.var.s.dotted$tissue)
  dat.var.s.dotted$tissue <- factor(as.character(dat.var.s.dotted$tissue), levels = c("Liver WT", "Liver KO", "Kidney WT", "Kidney KO"))
  dat.var.s.noisefloor.bytiss2 <- MakeTissueGeno(dat.var.s.noisefloor.bytiss)
  dat.var.s.noisefloor.bytiss2 <- UpdateStartEnd(dat.var.s.noisefloor.bytiss2)
  cbPalette.livkid <- c("#56B4E9", "grey90")
  
  # plot total variance
  plot.spectral.power.dotted <- ggplot() + 
    geom_bar(data=subset(dat.var.s.dotted, period.factor.cond != "Other"), mapping=aes(x = tissue, y = sum_sqr_mod, fill = period.factor.cond, linetype = geno),
             stat="identity", position=position_dodge(), colour="black", width = barwidth) +
    theme_bw(24) + theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1), 
                         panel.grid.major = element_blank(), 
                         panel.grid.minor = element_blank(), 
                         panel.background = element_blank(), 
                         axis.line = element_line(colour = "black"),
                         aspect.ratio = 0.5, 
                         legend.position="bottom") + xlab("") + ylab(jylab) +
    scale_y_continuous(breaks = c(seq(ymin, ymax, by = yinc)), limits = c(ymin, ymax)) + 
    scale_linetype_manual(values = c("solid", "dashed")) +
    scale_fill_manual(name = "Period [h]", drop=TRUE, values=cbPalette)
 
  plot.spectral.power.livkid <- ggplot() + 
    geom_bar(data=subset(dat.var.s2, period.factor == 24), mapping=aes(x = tissue, y = sum_sqr_mod, fill = geno),
             stat="identity", position=position_dodge(), colour="black", width = barwidth) +
    theme_bw(24) + theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1), 
                         panel.grid.major = element_blank(), 
                         panel.grid.minor = element_blank(), 
                         panel.background = element_blank(), 
                         axis.line = element_line(colour = "black"),
                         aspect.ratio = 0.5, 
                         legend.position="bottom") + xlab("") + ylab(jylab) +
    scale_y_continuous(breaks = c(seq(ymin, ymax, by = yinc)), limits = c(ymin, ymax)) + 
    scale_fill_manual(name = "Period [h]", drop=TRUE, values=cbPalette.livkid)
  
  
  # normalize to total variance
  dat.var.s2_adj <- dat.var.s.dotted %>%
    group_by(tissue) %>%
    mutate(s2_normalized = sum_sqr_mod / sum(sum_sqr_mod))
  
  dat.var.s2_adj.24 <- subset(dat.var.s2_adj, period.factor.cond == "24")
  dat.var.s2_adj.24$tissue <- factor(dat.var.s2_adj.24$tissue,
                                     levels = dat.var.s2_adj.24$tissue[order(dat.var.s2_adj.24$s2_normalized, decreasing = TRUE)])
  
  dat.var.s2_adj.12 <- subset(dat.var.s2_adj, period.factor.cond == "12")
  dat.var.s2_adj.12$tissue <- factor(dat.var.s2_adj.12$tissue,
                                     levels = dat.var.s2_adj.12$tissue[order(dat.var.s2_adj.12$s2_normalized, decreasing = TRUE)])
  dat.var.s2_adj$tissue <- factor(dat.var.s2_adj$tissue,
                                  levels = dat.var.s2_adj.24$tissue[order(dat.var.s2_adj.24$s2_normalized, decreasing = TRUE)])
  dat.var.s2_adj.noisefloor.bytiss <- subset(dat.var.s2_adj, period %in% noise.components) %>%
    group_by(tissue) %>%
    summarise(s2_normalized = mean(s2_normalized))
  dat.var.s2_adj.noisefloor.bytiss$start <- seq(lstart, by = 1, length.out = nrow(dat.var.s2_adj.noisefloor.bytiss))
  dat.var.s2_adj.noisefloor.bytiss$end <- seq(lstart + barwidth, by = 1, length.out = nrow(dat.var.s2_adj.noisefloor.bytiss))
  
  plot.normalized.spectral.power <- ggplot() + 
    geom_bar(data=subset(dat.var.s2_adj, period.factor.cond != "Other"), mapping=aes(x = tissue, y = s2_normalized, fill = period.factor.cond, linetype = geno),
             stat="identity", position=position_dodge(), colour="black", width = barwidth) +
    theme_bw(24) + theme(axis.text.x=element_text(angle=45,vjust = 1, hjust = 1), 
                         panel.grid.major = element_blank(), 
                         panel.grid.minor = element_blank(), 
                         panel.background = element_blank(), 
                         axis.line = element_line(colour = "black"),
                         aspect.ratio = 0.5, 
                         legend.position="bottom") + xlab("") + ylab("Fraction of Total Variance") +
    scale_linetype_manual(values = c("solid", "dashed")) +
    scale_fill_manual(name = "Period [h]", drop=TRUE, values=cbPalette)
  
  
  print(plot.spectral.power.dotted)
  print(plot.spectral.power.livkid)
  print(plot.spectral.power.livkid + geom_segment(data = dat.var.s.noisefloor.bytiss2, mapping=aes(x=start, xend=end, y=sum_sqr_mod.mean, yend=sum_sqr_mod.mean), linetype = jltype))
  print(plot.spectral.power.dotted + geom_segment(data = dat.var.s.noisefloor.bytiss, mapping=aes(x=start, xend=end, y=sum_sqr_mod.mean, yend=sum_sqr_mod.mean), linetype = jltype))
  print(plot.normalized.spectral.power <- plot.normalized.spectral.power + geom_segment(data = dat.var.s2_adj.noisefloor.bytiss, mapping=aes(x=start, xend=end, y=s2_normalized, yend = s2_normalized), linetype = jltype))
}


dev.off()