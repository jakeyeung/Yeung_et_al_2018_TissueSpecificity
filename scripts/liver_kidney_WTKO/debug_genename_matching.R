# 2016-07-12
# Jake Yeung
# why is matching
# Using 215/638 genes because they are in common
# and 
# "Using 157/638 genes because they are in common."
# "Using 300/1566 genes because they are in common."

# use same function to read mat to ram
source('/home/yeung/projects/ridge-regression/ridgeInR.R')


# Load --------------------------------------------------------------------

load("Robjs/liver_kidney_atger_nestle/fits.long.multimethod.filtbest.staggeredtimepts.Robj", v=T)

load("Robjs/S.long.multigene.filt.50000.Robj", v=T)
S.long <- subset(S.long, tissue %in% c("Kidney", "Liver"))
jcutoff <- 3
jcutoff.low <- 0

kidney.i <- 1
liver.i <- 2
S.sub.liverpeaks <- subset(S.long, gene == "Tars" & tissue %in% c("Kidney", "Liver")) %>%
  group_by(peak, gene) %>%  
  filter(zscore[liver.i] > jcutoff & zscore[kidney.i] < jcutoff.low)

exp <- "/home/yeung/projects/tissue-specificity/data/gene_exprs/liver_v_kidney/atger_with_kidney/Liver_SV129.mat"
site <- "/home/yeung/projects/tissue-specificity/data/sitecounts/motevo/liver_kidney_sitecounts_tissuespecpeaks/sitecounts_enhancers.model_Liver_SV129.method_g=1001.dist_40000.cutoff_3.cross_TRUE.tspeaks_TRUE.mat"

E <- read.table.handlerows(exp)
N <- read.table.handlerows(site)

# E = log2(as.matrix(E) + 1.0)  # convert to log2 BEFOREhand
E = as.matrix(E)
N = as.matrix(N)

# Es = center.cols(E)
# Ns = center.cols(N)

# take common rows
common.genes <- intersect(rownames(E), rownames(N))
max.genes <- max(nrow(E), nrow(N))
sprintf('Using %s/%s genes because they are in common.', length(common.genes), max.genes)
E <- E[common.genes, ]
N <- N[common.genes, ]



