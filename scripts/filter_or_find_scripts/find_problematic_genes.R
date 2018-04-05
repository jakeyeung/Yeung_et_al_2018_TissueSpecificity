# Jake Yeung
# 2015-09-18
# find_problematic_genes.R

MergeRows <- function(dat){
  # merge two rows into one to plot biological replicates
  dat.out <- data.frame(exprs1 = dat$exprs[1], exprs2 = dat$exprs[2])
  return(dat.out)
}

AssignIndex <- function(x){
  return(seq(length(x)))
}

# Source ------------------------------------------------------------------

source("scripts/functions/VarianceFunctions.R")
source("scripts/functions/ShannonEntropy.R")


# Load --------------------------------------------------------------------

load("Robjs/dat.long.fixed_rik_genes.Robj")


# Example -----------------------------------------------------------------

jgene <- "0610007P14Rik"
jgene <- "Svs1"
dat.sub <- subset(dat.long, gene == jgene)

dat.sub$time.m <- factor(dat.sub$time %% 24)

ggplot(dat.sub, aes(x = time.m, y = exprs, group = experiment, fill = experiment, colour = experiment)) + geom_point() + geom_line() + facet_wrap(~tissue) + ggtitle(jgene)

# Variation between replicates --------------------------------------------

dat.sub.var <- dat.sub %>%
  group_by(gene, tissue, time.m, experiment) %>%
  summarise(Var = var(exprs))

ggplot(dat.sub.var, aes(x = time.m, y = Var, group = experiment, fill = experiment)) + geom_bar(stat = "identity", position=position_dodge()) + facet_wrap(~tissue) + ggtitle(jgene)
ggplot(subset(dat.sub.var, experiment == "rnaseq"), aes(x = time.m, y = Var, group = experiment, fill = experiment)) + geom_bar(stat = "identity", position=position_dodge()) + facet_wrap(~tissue) + ggtitle(jgene)


# Find noisy genes genome-wide --------------------------------------------

dat.long$time.m <- dat.long$time %% 24

dat.24var <- subset(dat.long, experiment == "rnaseq") %>%
  group_by(gene, tissue, time.m, experiment) %>%
  summarise(Var = var(exprs))

jgene <- "Myh1"
ggplot(subset(dat.24var, gene == jgene), aes(x = time.m, y = Var)) + geom_bar(stat = "identity") + facet_wrap(~tissue)


# Take max and rank by noisiest -------------------------------------------

dat.maxvar <- dat.24var %>%
  group_by(gene, tissue, experiment) %>%
  summarise(max.Var = median(Var))

dat.maxvar <- dat.maxvar %>%
  group_by(tissue) %>%
  arrange(desc(max.Var)) %>%
  mutate(i = AssignIndex(max.Var))

# Plot distribution in histograms -----------------------------------------

# ggplot(dat.maxvar, aes(x = max.Var)) + geom_histogram(binwidth = 0.1) + facet_wrap(~tissue) + scale_x_log10()
# ggplot(dat.maxvar, aes(x = max.Var)) + stat_ecdf() + facet_wrap(~tissue)
ggplot(dat.maxvar, aes(x = max.Var, colour = tissue, fill = tissue)) + stat_ecdf()


# By entropy --------------------------------------------------------------

dat.entropy <- dat.24var %>%
  group_by(gene, tissue, experiment) %>%
  summarise(H = ShannonEntropy(p = Var, normalize = TRUE))

ggplot(dat.entropy, aes(x = H)) + geom_histogram(binwidth = 0.1) + facet_wrap(~tissue)
ggplot(dat.entropy, aes(x = H)) + stat_ecdf() + facet_wrap(~tissue)


# Plot biological replicates ----------------------------------------------

# 8 minutes
# dat.reps <- subset(dat.long, experiment == "rnaseq") %>%
#   group_by(gene, time.m, tissue, experiment) %>%
#   do(MergeRows(.))
# head(dat.reps)
load("Robjs/dat.reps.Robj", verbose = T)

# Plot replicates across tissues ------------------------------------------

dat.reps$time.m <- factor(dat.reps$time.m)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
ggplot(dat.reps, aes(x = exprs1, y = exprs2, colour = time.m)) + geom_point(alpha = 0.5) + facet_wrap(~tissue) + scale_colour_manual(values=cbPalette)


# Plot top 5 for each tissue ----------------------------------------------

outdir <- "plots/noisy_genes"
dir.create(outdir)
pdf(file.path(outdir, "noisy_genes.pdf"))

N <- 24
for (tiss in unique(dat.maxvar$tissue)){
  dat.sub <- subset(dat.maxvar, tissue == tiss & i <= N)
  genes.sub <- dat.sub$gene
  varmax <- signif(max(dat.sub$max.Var), digits = 3)
  varmin <- signif(min(dat.sub$max.Var), digits = 3)
  m <- ggplot(subset(dat.long, gene %in% genes.sub & tissue == tiss), aes(x = time, y = exprs, group = experiment, colour = experiment, fill = experiment)) + 
    geom_point() + geom_line() + facet_wrap(~gene) + ggtitle(paste(tiss, varmax, varmin))
  print(m)
}

dev.off()