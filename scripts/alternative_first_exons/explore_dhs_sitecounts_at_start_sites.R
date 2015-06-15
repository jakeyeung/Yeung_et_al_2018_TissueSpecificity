# Jake Yeung
# 2015-06-04
# explore_dhs_sitecounts_at_start_sites.R

rm(list=ls())

library(gplots)
library(reshape2)
library(dplyr)
library(ggplot2)
# Functions ---------------------------------------------------------------

source("scripts/functions/LoadSitecounts.R")

# motevo.long <- LoadSitecountsEncodeAll()
# save(motevo.long, file = "Robjs/motevo.long.Robj")
load("Robjs/motevo.long.Robj")
# Load outputs from filter_kallisto_for_start_sites.R
load(file = "Robjs/dat.rhyth.one.promoter.Robj")
load(file = "Robjs/tpm.afe.filt.Robj")
load(file = "Robjs/fit.afe.Robj")
load(file = "Robjs/dat.rhyth.Robj")
load(file = "Robjs/arrayrnaseq.mydat.Robj")



# Let's look at some known genes ------------------------------------------

jgene <- "Avpr1a"
jgene <- "Ube2u"
jgene <- "Insig2"
jensemblid <- "ENSMUST00000161068,ENSMUST00000071064"
jgene <- "Elovl3"
jgene <- "Osbpl8"
jgene <- "Tcap"
jgene <- "Asap2"
jgene <- "Tars"
jgene <- "Hn1l"
jgene <- "Rorc"

jgene <- "Tars"
motevo.sub <- subset(motevo.long, gene == jgene)

m <- dcast(data = motevo.sub, formula = motif ~ tissue, value.var = "motevo.value")
rownames(m) <- m$motif
m$motif <- NULL
m <- as.matrix(m)
# normalize by sum
m <- sweep(m, 2, colSums(m),"/")
m <- m[which(rowMeans(m) > 0), ]

par(mar=c(0, 12, 4.1, 2.1))
my_palette <- colorRampPalette(c("white", "black"))(n = 299)
heatmap.2(m, density.info = "density", trace = "none", margins = c(5, 10), main = jgene, col = my_palette, cexRow=0.4)

m.liver <- m[, "Liver"]
m.other <- m[, c("Cerebellum", "Heart", "Kidney", "Lung", "SkeletalMuscle")]
m.other.avg <- rowMeans(m.other)
abs.diff <- m.liver - m.other.avg
abs.diff <- abs.diff[order(abs.diff, decreasing = TRUE)]
(head(sort(abs.diff, decreasing = TRUE), n = 100))

par(mar=c(5, 12, 4.1, 2.1))
barplot(abs.diff[1:50], names.arg = names(abs.diff[1:50]), las = 2)

# Can we find enrichment of core clock motifs in tissue-specific r --------

# Normalize motevo.long
motevo.norm <- motevo.long %>%
  group_by(gene, ensemblid, tissue) %>%
  mutate(motevo.norm = motevo.value / sum(motevo.value))

dat.rhyth.one.promoter.summary <- dat.rhyth.one.promoter %>%
  group_by(gene) %>%
  summarise(min.pval=min(pval), max.amp=max(amp)) %>%
  arrange(min.pval)

top.n <- 200
top.hits <- dat.rhyth.one.promoter.summary$gene[1:top.n]

dat.rhyth.sub <- subset(dat.rhyth.one.promoter, gene %in% top.hits)
rhythmic.dic <- setNames(object = dat.rhyth.sub$is.rhythmic, 
                         nm = paste(dat.rhyth.sub$gene, dat.rhyth.sub$tissue, sep = ";"))

# First find enrichment of RORA.p2 as a test
# jmotifs <- c("RORA.p2", "NFIL3.p2", "bHLH_family.p2", "ATF2.p2", "SREBF1.2.p2", "HSF1.2.p2", "SRF.p3", "ZNF423.p2")
# jmotifs <- c("TFAP4.p2")
# jmotifs <- c("ESRRA.p2")
# jmotifs <- c("NFIL3.p2")

jmotifs <- sample(x = unique(motevo.long$motif), size = length(jmotifs), replace = FALSE)  # run as control
# TFAP4.p2 is interesting
# jmotifs.control <- c("ESRRA.p2")
# jmotifs.control <- c("HNF4A_NR2F1.2.p2")
# jmotifs.control <- c("ONECUT1.2.p2")

motevo.sub <- subset(motevo.norm, gene %in% top.hits & motif %in% jmotifs)

# motevo.sub.sum <- subset(motevo.norm, motif %in% jmotifs)
motevo.sub.sum <- motevo.sub %>%
  group_by(gene, tissue, ensemblid) %>%
  summarise(motevo.value.agg = sum(motevo.norm))

# Annotate motevo sitecounts to liver based on rhythnmicity
motevo.sub.sum$is.rhythmic <- rhythmic.dic[paste(motevo.sub.sum$gene, motevo.sub.sum$tissue, sep = ";")]

ggplot(data = motevo.sub.sum, aes(x = motevo.value.agg, fill = is.rhythmic)) + geom_density(alpha = 0.3) + ggtitle(paste(jmotifs, collapse = ", "))
ggplot(data = motevo.sub.sum, aes(x = is.rhythmic, y = motevo.value.agg)) + geom_boxplot()

# Are the two distributions different between TRUE and FALSE? -------------

x.true <- subset(motevo.sub.sum, is.rhythmic == TRUE)$motevo.value.agg
x.false <- subset(motevo.sub.sum, is.rhythmic == FALSE)$motevo.value.agg
ks.out <- ks.test(x.true, x.false)
print(ks.out)


# Loop for all motifs -----------------------------------------------------


# Do for all motifs
all.motifs <- unique(motevo.long$motif)
motevo.sub <- subset(motevo.norm, motif %in% all.motifs)
motevo.sub$is.rhythmic <- rhythmic.dic[paste(motevo.sub$gene, motevo.sub$tissue, sep = ";")]

# ks.outs <- motevo.norm %>%
#   filter(motif %in% all.motifs) %>%
#   group_by(motif) %>%
#   do(RunKsTest(.))

motevo.norm.sub <- subset(motevo.norm, motif %in% all.motifs)
ks.outs <- matrix(nrow = length(all.motifs), ncol = 2, dimnames = list(all.motifs, c("D", "ks.pval")))
for (jmotifs in all.motifs){
  motevo.sub <- subset(motevo.norm.sub, motif == jmotifs)
  motevo.sub$is.rhythmic <- rhythmic.dic[paste(motevo.sub$gene, motevo.sub$tissue, sep = ";")]
  x.true <- subset(motevo.sub, is.rhythmic == TRUE)$motevo.norm
  x.false <- subset(motevo.sub, is.rhythmic == FALSE)$motevo.norm
  ks.out <- ks.test(x.true, x.false)
  ks.outs[jmotifs, 1] <- ks.out$statistic
  ks.outs[jmotifs, 2] <- ks.out$p.value
}
ks.outs <- ks.outs[order(ks.outs[, 2]), ]
