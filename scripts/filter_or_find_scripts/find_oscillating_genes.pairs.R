# 2015-07-08 
# Find pairs between liver and other tissues. I hypothesize
# there will be more liver + kidney than liver + lung for example

library(dplyr)
# Functions ---------------------------------------------------------------

source("scripts/functions/LoadArrayRnaSeq.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/TissueSpecificRhythmicsFunctions.R")


# Load --------------------------------------------------------------------

dat.long <- LoadArrayRnaSeq()

load("Robjs/dat.rhyth.relamp.pvalmin1e-5.pvalmax0.05.relampmax.0.1.meancutoff6.Robj", verbose = T)

# how many rhythmic genes?
n.rhyth <- dat.rhyth.relamp %>%
  group_by(tissue) %>%
  summarise(n.rhyth = length(is.rhythmic[which(is.rhythmic == TRUE)]))


# How many are rhythmic specifically in liver + another tiss? -------------

n.pairs <- dat.rhyth.relamp %>%
  group_by(gene) %>%
  do(RhythTiss(.))

n.pairs.counts <- n.pairs %>%
  group_by(n.tiss) %>%
  summarise(count = length(n.tiss))

n.pairs.tissuecounts <- n.pairs %>%
  filter(n.tiss > 1) %>%
  group_by(tissues) %>%
  summarise(count = length(tissues)) %>%
  mutate(n.tiss = sapply(tissues, function(x) length(strsplit(x, split = ",")[[1]]))) %>%
  arrange(desc(count))

# global tissue-specific rhythmic
ggplot(subset(n.pairs.counts, n.tiss >= 1), aes(x = n.tiss, y = count)) + 
  geom_bar(stat = "identity") +
  scale_x_discrete() + scale_y_log10() + 
  ggtitle("Tissue-specific rhythmic genes: how many are rhythmic across different tissues") +
  ylab("Number of genes") +
  xlab("Number of tissues that is rhythmic for gene")


# How many rhythmic genes across tissues? ---------------------------------

n.rhyth$tissue <- factor(n.rhyth$tissue, 
                         levels = n.rhyth[order(n.rhyth$n.rhyth, decreasing = TRUE), ]$tissue)
ggplot(n.rhyth, aes(x = tissue, y = n.rhyth)) + geom_bar(stat = "identity")


# Let's look at pairs liver with others -----------------------------------

# pairs of genes: how are they? subset on liver
n.pairs.tiss.sub <- subset(n.pairs.tissuecounts, n.tiss == 2)
liver.pairs <- n.pairs.tiss.sub$tissues[grep("Liver", n.pairs.tiss.sub$tissues)]
n.pairs.liver <- subset(n.pairs.tiss.sub, tissues %in% liver.pairs)
n.pairs.liver$tissues <- factor(n.pairs.liver$tissues, levels = n.pairs.liver$tissues)
# reorder factors
n.pairs.tiss.sub$tissues <- factor(n.pairs.tiss.sub$tissues, levels = n.pairs.tiss.sub$tissues)

ggplot(n.pairs.liver, aes(x = factor(tissues), y = count)) + geom_bar(stat = "identity")


# How many Liver+Kidney+_ vs Liver+BFAT+_ ---------------------------------

liver.kid <- n.pairs.tissuecounts$tissues[grep("Kidney.*Liver", n.pairs.tissuecounts$tissues)]
liver.bfat <- n.pairs.tissuecounts$tissues[grep("BFAT.*Liver", n.pairs.tissuecounts$tissues)]
n.liverkid <- subset(n.pairs.tissuecounts, tissues %in% liver.kid)
n.liverbfat <- subset(n.pairs.tissuecounts, tissues %in% liver.bfat)
# sum them up
n.sum.liverkid <- sum(n.liverkid$count)
n.sum.liverbfat <- sum(n.liverbfat$count)


# Get kidney & liver genes ------------------------------------------------

# TODO: Make this into a function yo

genes.kidliv <- subset(n.pairs, tissues == "Kidney,Liver")
genes.bfatliv <- subset(n.pairs, tissues == "BFAT,Liver")
genes.aortabfat <- subset(n.pairs, tissues == "Aorta,BFAT")
genes.bfat <- subset(n.pairs, tissues == "BFAT")
# genes.bfatmus <- subset(n.pairs, tissues == "BFAT,Mus")

phases.kidliv <- subset(dat.rhyth.relamp, gene %in% genes.kidliv$gene & tissue %in% c("Kidney", "Liver"))
phases.bfatliv <- subset(dat.rhyth.relamp, gene %in% genes.bfatliv$gene & tissue %in% c("BFAT", "Liver"))
phases.aortabfat <- subset(dat.rhyth.relamp, gene %in% genes.aortabfat$gene & tissue %in% c("Aorta", "BFAT"))
phases.bfat <- subset(dat.rhyth.relamp, gene %in% genes.bfat$gene & tissue %in% c("BFAT"))

phases.kidliv.diff <- phases.kidliv %>%
  group_by(gene) %>%
  summarise(phasediff = GetPhaseDiff(phase))

phases.bfatliv.diff <- phases.bfatliv %>%
  group_by(gene) %>%
  summarise(phasediff = GetPhaseDiff(phase))

phases.aortabfat.diff <- phases.aortabfat %>%
  group_by(gene) %>%
  summarise(phasediff = GetPhaseDiff(phase))  # liver minus bfat

ggplot(phases.kidliv.diff, aes(x = phasediff)) + 
  geom_histogram(binwidth = 1) + ggtitle("Liver phase - Kidney phase")  # liver minus kidney shows coherence

ggplot(phases.bfatliv.diff, aes(x = phasediff)) + 
  geom_histogram(binwidth = 1) + ggtitle("Liver phase - BFAT phase")  # liver minus bfat show discordance

ggplot(phases.aortabfat.diff, aes(x = phasediff)) + 
  geom_histogram(binwidth = 1)

ggplot(phases.bfat, aes(x = phase)) + geom_histogram(binwidth = 0.5) + coord_polar(theta = "x") + scale_x_continuous(limits = c(0, 24), breaks = seq(2, 24, 2))

ggplot(phases.bfat, aes(x = phase, y = relamp)) + geom_point()

# Phase distribution of liver genes vs brown fat rhythmic genes -----------

ggplot(subset(phases.bfatliv), aes(x = phase, colour = tissue, fill = tissue)) + 
  geom_histogram(alpha = 0.4, binwidth = 0.5) + coord_polar(theta = "x") +
  scale_x_continuous(limits = c(0, 24), breaks = seq(2, 24, 2))

ggplot(subset(phases.kidliv), aes(x = phase, colour = tissue, fill = tissue)) + 
  geom_histogram(alpha = 0.4, binwidth = 0.5) + coord_polar(theta = "x") +
  scale_x_continuous(limits = c(0, 24), breaks = seq(2, 24, 2))

ggplot(phases.aortabfat, aes(x = phase, colour = tissue, fill = tissue)) + 
  geom_histogram(alpha = 0.4, binwidth = 0.5) + coord_polar(theta = "x") +
  scale_x_continuous(limits = c(0, 24), breaks = seq(2, 24, 2))


# Phase distribution of rhythmic genes: Liver vs BFAT ---------------------

rhythmic.genes <- dat.rhyth.relamp %>%
  group_by(tissue) %>%
  do(GetRhythmic(.))

# plot phase of all rhythmic genes
ggplot(rhythmic.genes, aes(x = phase)) +
  geom_histogram(alpha = 0.8, binwidth = 0.5) + coord_polar(theta = "x") +
  scale_x_continuous(limits = c(0, 24), breaks = seq(2, 24, 2)) +
  facet_wrap(~tissue)

# plot phase without liver
ggplot(subset(rhythmic.genes, tissue == "BFAT"), aes(x = phase)) +
  geom_histogram(alpha = 0.8, binwidth = 0.5) + coord_polar(theta = "x") +
  scale_x_continuous(limits = c(0, 24), breaks = seq(2, 24, 2)) +
  facet_wrap(~tissue)


# Phase distribution of tissue-specific rhythmic genes --------------------

tissue.spec.rhyth <- subset(n.pairs, n.tiss == 1)

tissue.spec.rhyth.genes <- dat.rhyth.relamp %>%
  group_by(tissue) %>%
  do(GetTissueSpec(., tissue.spec.rhyth))

ggplot(subset(tissue.spec.rhyth.genes, tissue != "Liver"), aes(x = phase)) +
  geom_histogram(alpha = 0.8, binwidth = 0.5) + coord_polar(theta = "x") +
  scale_x_continuous(limits = c(0, 24), breaks = seq(2, 24, 2)) +
  facet_wrap(~tissue)