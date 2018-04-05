# 2015-09-24
# Jake Yeung
# plot_core_clock_genes_disruption.R
# are core clock genes differing between different tissues


# Function ----------------------------------------------------------------

source("scripts/functions/GetClockGenes.R")
clockgenes <- GetClockGenes()


# Load long dat -----------------------------------------------------------

load("Robjs/dat.long.fixed_rik_genes.Robj", verbose = T)


# Fit ---------------------------------------------------------------------

dat.sub <- subset(dat.long, gene %in% clockgenes)

per <- 24
w <- 2 * pi / per
dat.fit <- dat.sub %>%
  group_by(gene, tissue) %>%
  do(mod = lm(exprs ~ 0 + experiment + cos(w * time) + sin(w * time), data = .))

dat.fit$var <- apply(dat.fit, MARGIN = 1, function(x) summary(x[3][[1]])$coefficients[["cos(w * time)", "Std. Error"]] ^ 2)

dat.fit$cos <- apply(dat.fit, MARGIN = 1, function(x) coef(x[3][[1]])[["cos(w * time)"]])
dat.fit$coslow <- apply(dat.fit, MARGIN = 1, function(x) confint(x[3][[1]])[[3, 1]])
dat.fit$coshigh <- apply(dat.fit, MARGIN = 1, function(x) confint(x[3][[1]])[[3, 2]])
dat.fit$sin <- apply(dat.fit, MARGIN = 1, function(x) coef(x[3][[1]])[["sin(w * time)"]])
dat.fit$sinlow <- apply(dat.fit, MARGIN = 1, function(x) confint(x[3][[1]])[[4, 1]])
dat.fit$sinhigh <- apply(dat.fit, MARGIN = 1, function(x) confint(x[3][[1]])[[4, 2]])

ggplot(dat.fit, aes(x = cos, y = sin, colour = gene, xmin = coslow, xmax = coshigh, ymin = sinlow, ymax = sinhigh)) + geom_point() + 
  scale_x_continuous(limits = c(-5, 5)) + 
  scale_y_continuous(limits = c(-5, 5)) + 
  facet_wrap(~tissue) + 
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_errorbar() + geom_errorbarh()

ggplot(subset(dat.fit, gene %in% c("Arntl", "Nr1d1", "Cry1")), aes(x = cos, y = sin, label = tissue, colour = gene, xmin = coslow, xmax = coshigh, ymin = sinlow, ymax = sinhigh)) + 
  geom_point(size = 5) + 
  scale_x_continuous(limits = c(-5, 5)) + 
  scale_y_continuous(limits = c(-5, 5)) + 
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_errorbar() + geom_errorbarh() + 
  geom_text(colour = "black")

ggplot(subset(dat.fit, gene %in% c("Arntl", "Nr1d1", "Cry1")), aes(x = cos, y = sin, label = tissue, colour = gene, xmin = coslow, xmax = coshigh, ymin = sinlow, ymax = sinhigh)) + 
  geom_point(size = 5) + 
  scale_x_continuous(limits = c(-5, 5)) + 
  scale_y_continuous(limits = c(-5, 5)) + 
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_errorbar() + geom_errorbarh() + 
  geom_text(colour = "black") +
  facet_wrap(~tissue)


