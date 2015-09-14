# Jake Yeung
# 2015-09-14
# standardize_to_reverba.R
# Perform SVD after Fourier to reverba

ref.gene <- "Nr1d1"

# Functions ---------------------------------------------------------------

library(hash)
source("scripts/functions/SvdFunctions.R")
source("scripts/functions/LoadArrayRnaSeq.R")

# Load data from heatmap_tissue_specific_rhythms.R ------------------------

dat.long <- LoadArrayRnaSeq(fix.rik.xgene = TRUE)
load(file = "Robjs/dat.fit.Robj", verbose = T)


# Get relamp --------------------------------------------------------------

dat.fit.relamp <- GetRelampByGene(dat.fit, by.gene = ref.gene)

dat.fit.relamp <- dat.fit.relamp %>%
  group_by(gene) %>%
  mutate(max.exprs = max(int.rnaseq), min.pval = min(pval), max.relamp = max(relamp), relamp.norm = relamp / sum(relamp), amp.norm = amp / sum(amp))

# Get dat.complex ---------------------------------------------------------

genes.exprs <- unique(subset(dat.fit.relamp, max.exprs >= 4)$gene)

dat.complex <- TemporalToFrequencyDatLong(subset(dat.long, gene %in% genes.exprs), period = 24, n = 8, interval = 6, add.entropy.method = "array")

# z-score
dat.complex$z <- dat.complex$exprs.transformed * dat.complex$frac.weight
dat.complex$mod.z <- Mod(dat.complex$z)

# Adjust to nr1d1 but not to z-score
ref.amps <- subset(dat.complex, gene == ref.gene, select = c(tissue, gene, exprs.transformed))
ref.amps.dic <- hash(ref.amps$tissue, Mod(ref.amps$exprs.transformed))

dat.complex$exprs.adj <- mapply(function(tiss, exprs) exprs / ref.amps.dic[[tiss]], as.character(dat.complex$tissue), dat.complex$exprs.transformed)
dat.complex$mod.exprs.adj <- Mod(dat.complex$exprs.adj)

# adjust to nr1d1 by z-score
ref.amps.z <- subset(dat.complex, gene == ref.gene, select = c(tissue, gene, z))
ref.amps.z.dic <- hash(ref.amps$tissue, Mod(ref.amps.z$z))

dat.complex$exprs.adj.z <- mapply(function(tiss, exprs) exprs / ref.amps.z.dic[[tiss]], as.character(dat.complex$tissue), dat.complex$z)
dat.complex$mod.exprs.adj.z <- Mod(dat.complex$exprs.adj.z)

# SVD on all genes --------------------------------------------------------

# filt.tiss <- c("Cere", "Hypo", "BS")
filt.tiss <- c()
s.all <- SvdOnComplex(subset(dat.complex, ! tissue %in% filt.tiss), value.var = "exprs.transformed")
s.all.z <- SvdOnComplex(subset(dat.complex, ! tissue %in% filt.tiss), value.var = "z")
s.all.nr1d1 <- SvdOnComplex(subset(dat.complex, ! tissue %in% filt.tiss), value.var = "exprs.adj")
s.all.nr1d1.z <- SvdOnComplex(subset(dat.complex, ! tissue %in% filt.tiss), value.var = "exprs.adj.z")


# Plot --------------------------------------------------------------------

comp <- 1

eigens.all <- GetEigens(s.all, period = 24, comp = comp, adj.mag = TRUE)
eigens.all.z <- GetEigens(s.all.z, period = 24, comp = comp, adj.mag = TRUE)
eigens.all.nr1d1 <- GetEigens(s.all.nr1d1, period = 24, comp = comp, adj.mag = TRUE)
eigens.all.nr1d1.z <- GetEigens(s.all.nr1d1.z, period = 24, comp = comp, adj.mag = TRUE)

jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)

multiplot(eigens.all$v.plot, eigens.all$u.plot, layout = jlayout)
multiplot(eigens.all.z$v.plot, eigens.all.z$u.plot, layout = jlayout)
multiplot(eigens.all.nr1d1$v.plot, eigens.all.nr1d1$u.plot, layout = jlayout)
multiplot(eigens.all.nr1d1.z$v.plot, eigens.all.nr1d1.z$u.plot, layout = jlayout)


# Centering ---------------------------------------------------------------

# z-score is best
dat.complex <- dat.complex %>%
  group_by(gene) %>% 
  mutate(z.center = scale(z, center = TRUE, scale = FALSE),
         mod.z.center = Mod(z.center))
dat.complex$mod.exprs.adj <- NULL
dat.complex$mod.exprs.adj.z <- NULL

comp <- 3
s.all.z.center <- SvdOnComplex(subset(dat.complex, ! tissue %in% filt.tiss), value.var = "z.center")
eigens.all.z.center <- GetEigens(s.all.z.center, period = 24, comp = comp, adj.mag = TRUE)
multiplot(eigens.all.z.center$v.plot, eigens.all.z.center$u.plot, layout = jlayout)

