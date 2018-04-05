# 2015-07-10
# Look at Mef2c and other sitecounts in promoters only


# Functions ---------------------------------------------------------------

library(reshape2)
source("scripts/functions/LoadSitecounts.R")

# Explore sitecounts ------------------------------------------------------

N.list <- LoadSitecounts(gene_ids=FALSE)  # list of N and N.promoter
N.annot <- LoadEnsemblToPromoter()
# save(N.annot, file="Robjs/N.annot.Robj")
# rownames(N.annot) <- N.annot$ensemblid
N <- as.matrix(N.list$N)
N.promoter <- N.list$N.promoter

N.long <- LoadSitecountsPromotersLong()

# mef2c promoter IDs
mef2c_ids <- c("chr3:87945759-87946850(+)", "chr3:87945946-87947042(+)", "chr3:87973434-87974522(+)")
mef2c_ids <- c("chr13:83642532-83643532(+)", "chr13:83663317-83664504(+)", "chr13:83803330-83804349(+)")
N.sub <- subset(N.long, promoterid %in% mef2c_ids)

M <- dcast(data = N.sub, motif ~ promoterid, value.var = "sitecount")
rownames(M) <- M$motif
M$motif <- NULL

M <- t(scale(t(M), center = TRUE, scale = FALSE))

s <- prcomp(M)
screeplot(s)
biplot(s)

pca1 <- sort(s$x[, 1], decreasing = TRUE)
par(mar=c(15.1, 4.1, 4.1, 2.1))
barplot(pca1, cex.names = 0.7, las = 2)
