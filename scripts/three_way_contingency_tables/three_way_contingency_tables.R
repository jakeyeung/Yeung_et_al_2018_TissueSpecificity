# 2016-10-06
# Jake Yeung
# 3-way contingency tables


# Load --------------------------------------------------------------------

load("Robjs/three_way_cooccurence/three.way.cooccurrence.Robj", verbose = T)

fits.long <- do.call(rbind, fits)

fits.long <- fits.long[order(fits.long$pval), ]

# get pairs of RORA
tf <- "RORA.p2"
tf <- "bHLH_family.p2"
tf <- "ONECUT1.2.p2"
tf <- "FOXA2.p3"
tf <- "SRF.p3"
tf <- "HSF1,2.p2"
head(fits.long[grepl(tf, fits.long$pair), ], n = 50)
