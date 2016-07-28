# 2016-07-28
# See if we can find enrichment of motifs in promoters between kidney rhythmic and flat genes

rm(list=ls())

library(dplyr)
library(ggplot2)
source("scripts/functions/LdaFunctions.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/GetTFs.R")
source("scripts/functions/PlotFunctions.R")
source("scripts/functions/LiverKidneyFunctions.R")



# Regulation with LDA -----------------------------------------------------

# from LDA directory: sitecounts_analysis_promoter.R

load("Robjs/N.long.promoters_500.Robj", verbose=T)
# load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.15.phase_sd_maxdiff_avg.Robj", verbose=T)
load("Robjs/liver_kidney_atger_nestle/fits.long.multimethod.filtbest.staggeredtimepts.bugfixed.Robj", v=T)
load("Robjs/liver_kidney_atger_nestle/dat.long.liverkidneyWTKO.bugfixed.Robj", v=T)

dat.orig <- dat.long
dat.long <- CollapseTissueGeno(dat.long)
dat.long <- StaggeredTimepointsLivKid(dat.long)

dat.mean <- dat.long %>%
  group_by(gene, tissue, geno) %>%
  summarise(exprs.mean = mean(exprs))

fits.long.filt <- subset(fits.long.filt, method == "g=1001")


# Define gene sets --------------------------------------------------------

fits.count <- fits.long.filt %>%
  group_by(model) %>%
  summarise(count = length(gene)) %>%
  arrange(desc(count))


fg.models <- "Kidney_SV129,Kidney_BmalKO"
fg.models <- "Kidney_SV129"
fg.models <- "Liver_SV129"
flat.models <- ""
tissuewide.models <- "Liver_SV129,Kidney_SV129,Liver_BmalKO,Kidney_BmalKO"
genes.flat <- as.character(subset(fits.long.filt, model == flat.models)$gene)
genes.kid <- as.character(subset(fits.long.filt, model == fg.models)$gene)
genes.tw <- as.character(subset(fits.long.filt, model == tissuewide.models)$gene)

# show tissue-specificity
jmods <- c(fg.models, tissuewide.models, flat.models)
jmeth <- "g=1001"
for (jmod in jmods){
  print(jmod)
  genes <- as.character(subset(fits.long.filt, model == jmod & method == jmeth)$gene)
  print(PlotMeanExprsOfModel(dat.mean, genes = genes, jmodel = jmod, sorted = TRUE, avg.method = "mean"))
}


fg.mat <- LongToMat.lda(subset(N.long, gene %in% genes.kid))
flat.mat <- LongToMat.lda(subset(N.long, gene %in% genes.flat))
rhyth.mat <- LongToMat.lda(subset(N.long, gene %in% genes.tw))

out.flat <- DoLdaPromoters(fg.mat, flat.mat)
out.rhyth <- DoLdaPromoters(fg.mat, rhyth.mat)

out.df <- data.frame(motif = colnames(out.flat$x), discrim.flat = out.flat$discrim, discrim.rhyth = out.rhyth$discrim)

m <- ggplot(out.df, aes(x = discrim.flat, y = discrim.rhyth, label = motif)) + geom_text() + ggtitle(paste("Rhyth model:", paste(fg.models, collapse='|'))) + 
  geom_hline(aes(yintercept=0)) + geom_vline(aes(xintercept=0)) +
  theme_bw(24)
print(m)
PlotLdaOut(out.flat, jtitle = paste("Rhyth model:", paste(fg.models, collapse='|')), jcex = 0.5)
PlotSeparation(out.flat, jtitle = paste("Rhyth model:", paste(fg.models, collapse='|')))

# Show motifs of lowest projection (most liver like)
xproj <- out.flat$xproj
rownames(xproj) <- rownames(out.flat$x)

head(xproj[order(xproj), ])


# Plot expression of genes matching to motif  -----------------------------

dis <- out.flat$discrim
names(dis) <- colnames(out.flat$x)
dis <- dis[order(dis[, 1])]
# take top n
top.n <- 10
motifs <- names(head(dis, n = top.n))
tfs <- GetTFs(get.mat.only = TRUE)

for (jmotif in motifs){
  genes <- GetGenesFromMotifs(jmotif, tfs)
  for (g in genes){
    jsub <- subset(dat.long, gene == g)
    
    if (nrow(jsub) > 0){
      print(PlotGeneTissuesWTKO(jsub, jtitle = g))
    } else {
      print(paste("Could not find gene:", g))
    }
  }
}

