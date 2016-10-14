# 2016-09-20
# Jake Yeung
# Do motif enrichment on ~900 tissue wide genes (obtained from hogenesch) and
# ask whether these genes are enriched with HIC1, Ebox, Dbox, ROR, HSF, SRF, NR3C1.

rm(list=ls())

library(dplyr)
library(hash)
library(reshape2)
library(ggplot2)
library(ggrepel)
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/HandleMotifNames.R")
source("scripts/functions/SvdFunctions.R")
source("scripts/functions/LoadActivitiesLong.R")
source("scripts/functions/PlotActivitiesFunctions.R")
source("scripts/functions/CosSineFunctions.R")

# Load hogenesch hits -----------------------------------------------------

jmeth <- "g=1001"
load("Robjs/nconds_g1000_11_tissues/fits_long.11_tiss_3_max.g1000.bestmodel.filteramp.0.15.Robj", v=T)
load("Robjs/liver_kidney_atger_nestle/fits.long.multimethod.filtbest.staggeredtimepts.bugfixed.annotated.Robj", v=T)
fits.long.filt <- subset(fits.long.filt, method == jmeth)
load("Robjs/liver_kidney_atger_nestle/dat.long.liverkidneyWTKO.bugfixed.Robj", v=T); dat.wtko <- dat.long; rm(dat.long)
load("Robjs/liver_kidney_atger_nestle/fits.bytiss.bugfixed.Robj", v=T)
load("Robjs/N.long.promoters_500.Robj", v=T)
load("Robjs/dat.long.fixed_rik_genes.Robj", v=T)

# Explore -----------------------------------------------------------------

genes.tw <- as.character(subset(fits.long, n.rhyth >= 8)$gene)

# load MARA motifs

# get regulators: hogenesch 
omega <- 2 * pi / 24
zscore.min <- 1.25
N <- 4
outbase <- "/home/yeung/projects/tissue-specificity/results/MARA.hogenesch"
outmain <- file.path(outbase, paste0("promoters.tissuewide.filteramp.0.15.mat"))
indir <- file.path(outmain, "expressed_genes_deseq_int.centeredTRUE")
act.long <- LoadActivitiesLong(indir, shorten.motif.name = TRUE)
act.complex <- ProjectWithZscore(act.long, omega, n = N)
sig.motifs <- unique(as.character(subset(act.complex, zscore > zscore.min)$gene))
s.act <- SvdOnComplex(subset(act.complex, tissue != "WFAT" & gene %in% sig.motifs), value.var = "exprs.transformed")

comp <- 1; dotsize <- 2
eigens.act <- GetEigens(s.act, period = 24, comp = comp, adj.mag = TRUE, constant.amp = dotsize, label.n = 25, jtitle = "", peak.to.trough = TRUE, dot.col = "black", dotsize = 2, dotshape = 18)
print(eigens.act$u.plot)
print(eigens.act$v.plot)
multiplot(eigens.act$u.plot, eigens.act$v.plot, cols = 2)

# Filter out flat genes in Gachon -----------------------------------------

amp.cutoff <- 0
fits.sub <- subset(fits.long.filt, gene %in% genes.tw & model != "" & amp.avg > amp.cutoff) 

fits.sum <- fits.sub %>% 
  group_by(model) %>% 
  summarise(count = length(gene)) %>%
  arrange(desc(count))


# Identify clock and system driven models ---------------------------------

clock.or.system <- sapply(fits.sum$model, function(m){
  # if model is Liver_SV129,Kidney_SV129;Liver_BmalKO,Kidney_BmalKO, consider
  # it clock
  if (m == "Liver_SV129,Kidney_SV129;Liver_BmalKO,Kidney_BmalKO"){
    return("clock")
  }
  
  if (grepl("KO", m)){
    return("system")
  } else if (grepl("SV129", m)){
    return("clock")
  } else {
    warning("Neither KO or SV129")
    return(NA)
  }
})
clock.sys.hash <- hash(as.character(fits.sum$model), clock.or.system)

fits.sub$clksys <- sapply(as.character(fits.sub$model), function(m) clock.sys.hash[[m]])
clock.sys.gene.hash <- hash(as.character(fits.sub$gene), fits.sub$clksys)

# Penalized LDA between clock and system driven ----------------------------

motifs.tmp <- as.character(unique(N.long$motif))
motifs.hash <- hash(motifs.tmp, sapply(motifs.tmp, RemoveP2Name))

N.long$motif2 <- sapply(as.character(N.long$motif), function(m) motifs.hash[[m]])

# motifs <- c("HIC1", "RORA", "SRF", "NR3C1", "bHLH_family", "ATF2", "NFIL3", "HSF1.2",
            # "TFDP1", "NRF1", "SRY", "IRF1.2.7", "EP300", "SPZ1", "NFE2L2", "NFIX", "TEAD1")
# motifs <- unique(N.long$motif2)
motifs <- sig.motifs
# motifs <- c("HIC1", "RORA", "SRF", "NR3C1", "bHLH_family", "NFIL3", "HSF1.2")
# motifs <- c("HIC1", "RORA", "SRF", "NR3C1", "bHLH_family", "NFIL3", "HSF1.2", "ATF2", "NFIX", "TEAD1", "NFE2L2")
genes <- as.character(fits.sub$gene)
clksys <- as.character(fits.sub$clksys)

N.sub <- subset(N.long, motif2 %in% motifs & gene %in% genes)
N.sub$clksys <- sapply(as.character(N.sub$gene), function(g) clock.sys.gene.hash[[g]])

M.full <- dcast(N.sub, formula = gene + clksys ~ motif2, value.var = "sitecount", fun.aggregate = sum)

M <- subset(M.full, select = c(-gene, -clksys)); rownames(M) <- M.full$gene
M.labs <- as.numeric(as.factor(M.full$clksys))
# Do penalized LDA --------------------------------------------------------

library(penalizedLDA)

# remove columns that are all 0
M <- M[, which(colSums(M) != 0)]

jlambda <- 0.035  # liv only
jlambda <- 0.055  # liv only
plda.out <- PenalizedLDA(M, M.labs, lambda = jlambda, K = 1, standardized = FALSE)

# plot pretty
vec.length <- sqrt(plda.out$discrim[, 1]^2)

jsize.cutoff <- 0
jsize.pairs.cut <- sapply(vec.length, function(jsize){
  if (jsize > jsize.cutoff){
    return(jsize)
  } else {
    return(0)
  }
})

labels <- names(plda.out$x)
labels.cut <- mapply(function(jlab, jsize){
  if (jsize <= 0){
    return("")
  } else {
    return(jlab)
  }
}, labels, jsize.pairs.cut)

dat.plot <- data.frame(discrim = plda.out$discrim[, 1],
                       motif = labels.cut,
                       motif.orig = labels,
                       vec.length = vec.length,
                       vec.length.cut = jsize.pairs.cut) %>% 
  mutate(discrim.floor = ifelse(discrim > 0, "Systemic", "Clock")) %>%
  arrange(discrim) %>%
  mutate(Index = seq(length(discrim)))
  
dat.labs <- subset(dat.plot, vec.length.cut > 0)

m <- ggplot(dat.plot, aes(x = discrim.floor, y = discrim, label = motif)) + 
  geom_point(size = 0.01) + 
  geom_text_repel(size = 5) + 
  theme_bw() + 
  theme(aspect.ratio = 0.33, legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  xlab("") + ylab("Motif loadings") + 
  theme(aspect.ratio = 1)
print(m)

m.index <- ggplot(dat.plot, aes(x = Index, y = discrim, label = motif)) + 
  geom_point(size = 0.01) + 
  geom_text_repel(size = 5) + 
  theme_bw() + 
  theme(aspect.ratio = 0.33, legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  xlab("") + ylab("Motif loadings") + 
  theme(aspect.ratio = 1, 
        axis.ticks.x = element_blank())
print(m.index)

gene.plot <- data.frame(proj = plda.out$xproj, 
                        gene = rownames(plda.out$x),
                        jlabel = ifelse(plda.out$y == 1, "Clock", "Systems"))

mm <- ggplot(gene.plot, aes(y = proj, x = as.factor(jlabel), label = gene)) + 
  geom_boxplot() +
  geom_text() + 
  theme_bw() + 
  xlab("") + 
  ylab("Projection") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(mm)

# label colours by systems-driven, clock-driven, or "gray"

d.pos <- "red"
d.neg <- "blue"
d.zero <- "gray85"

s.pos <- "8"
s.neg <- "18"
s.zero <- "46"
col.hash <- hash(dat.plot$motif.orig, sapply(dat.plot$discrim, function(d){
  if (d > 0) return(d.pos)
  if (d < 0) return(d.neg)
  if (d == 0) return(d.zero)
}))

eigens.act <- GetEigens(s.act, period = 24, comp = comp, adj.mag = TRUE, 
                        constant.amp = dotsize, 
                        label.n = Inf, jtitle = "", 
                        peak.to.trough = TRUE, 
                        # dot.col = "black", 
                        dot.col = col.hash, 
                        dotsize = 4, 
                        dotshape = 18,
                        disable.text = TRUE)
print(eigens.act$u.plot)

# Do enrichment on a single motif -----------------------------------------

jmotif <- "NR3C1"

jmotifs <- motifs
for (jmotif in jmotifs){
  M.motif <- subset(M.full, select = c("gene", "clksys", jmotif))
  colnames(M.motif) <- c("gene", "clksys", "sitecount")
  print(ggplot(M.motif, aes(x = clksys, y = sitecount, label = gene)) + geom_boxplot() + geom_text() + ggtitle(jmotif) + theme_bw())
  print(paste("T test for motif:", jmotif))
  print(t.test(sitecount ~ clksys, M.motif))
}


# Save for reproducibility  -----------------------------------------------

save(M, M.labs, plda.out, file = "Robjs/liver_kidney_atger_nestle/systems_clockdriven_tissuewide_genes.Robj")

