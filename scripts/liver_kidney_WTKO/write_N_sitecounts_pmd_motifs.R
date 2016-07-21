# 2016-07-20
# Jake Yeung
# write_N_sitecounts_pmd_motifs.R
# Many cross products are correlated. Run penalized matrix decomposition to find
# loadings to map 190 motifs onto a 2 dimensional motif vector

rm(list=ls())

descrip <- "sep_liv_rhyth.RORA_bHLH_SRF.fixsignbugfixed"
K <- 3

setwd("~/projects/tissue-specificity")

library(ggplot2)
library(ggrepel)
library(hash)
library(reshape2)
library(hash)
library(PMA)


source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/PlotFunctions.R")
source("scripts/functions/SitecountsFunctions.R")
source("scripts/functions/LdaFunctions.R")
source("scripts/functions/RemoveP2Name.R")
source("scripts/functions/GetTopMotifs.R")
source("scripts/functions/HandleMotifNames.R")
source("scripts/functions/ListFunctions.R")
source("scripts/functions/PlotActivitiesFunctions.R")
source("scripts/functions/LoadActivitiesLong.R")
source("scripts/functions/LiverKidneyFunctions.R")
source("scripts/functions/SvdFunctions.R")
source("scripts/functions/NcondsAnalysisFunctions.R")
source("scripts/functions/FisherTestSitecounts.R")
source("scripts/functions/ColorFunctions.R")

PcToMotif <- function(v, pc, top.n = 3){
  x <- v[, pc]
  x <- x[which(x != 0)]
  motif.names <- names(sort(x, decreasing = TRUE))
  return(paste(motif.names[1:top.n], collapse = "-"))
}


# Load --------------------------------------------------------------------

pldarobj <- "/home/yeung/projects/tissue-specificity/Robjs/liver_kidney_atger_nestle/plda_robjs/2D.posterior.multigene.distfilt.morenonliv.bugfixed.liverWTKO.40000.cutoff.3.cutofflow0.method.g=1001.Robj"
# pldarobj <- "/home/yeung/projects/tissue-specificity/Robjs/liver_kidney_atger_nestle/plda_robjs/2D.posterior.multigene.distfilt.morenonliv.bugfixed.liverWTKO.40000.cutoff.3.cutofflow0.5.method.g=1001.Robj"
load(pldarobj)

load("Robjs/liver_kidney_atger_nestle/dat.long.liverkidneyWTKO.bugfixed.Robj")
dat.orig <- dat.long
dat.long <- CollapseTissueGeno(dat.long, keep.tissue.col = TRUE)
dat.long <- StaggeredTimepointsLivKid(dat.long)
load("Robjs/fits.relamp.Robj")
load("Robjs/liver_kidney_atger_nestle/fits.bytiss.Robj")
fits.bytiss <- subset(fits.bytiss, !is.na(gene))

jsub.hog <- subset(fits.relamp, tissue == "Liver")
phase.hog <- hash(as.character(jsub.hog$gene), jsub.hog$phase)
jsub.atg <- subset(fits.bytiss, tissue == "Liver_SV129") 
phase.atg <- hash(as.character(jsub.atg$gene), jsub.atg$phase)
amp.atg <- hash(as.character(jsub.atg$gene), jsub.atg$amp)
pval.atg <- hash(as.character(jsub.atg$gene), jsub.atg$pval)



# Setup mat ---------------------------------------------------------------

mat.fgbg.lab.lst.3 <- SetUpMatForLda(mat.fg, mat.bgnonliver, mat.bg, has.peaks = TRUE)
mat.fgbg.3 <- mat.fgbg.lab.lst.3$mat.fgbg; labels3 <- mat.fgbg.lab.lst.3$labels
colnames(mat.fgbg.3) <- sapply(colnames(mat.fgbg.3), RemoveP2Name)
colnames(mat.fgbg.3) <- sapply(colnames(mat.fgbg.3), function(cname){
  return(RemoveCommasBraces(cname))
}, USE.NAMES = FALSE)


rhyth.motifs <- sapply(GetTopMotifs("rhythmic"), RemoveP2Name, USE.NAMES = FALSE)
rhyth.motifs <- c(rhyth.motifs, c("SRF"))
# add DBP
# rhyth.motifs <- c(rhyth.motifs, c("DBP"))
rhyth.motifs <- rhyth.motifs[which( ! rhyth.motifs %in% c("ATF6"))]
tissue.motifs <- sapply(GetTopMotifs("tissue"), RemoveP2Name, USE.NAMES = FALSE)
tissue.motifs <- c(tissue.motifs, c("ATF5_CREB3", "ATF6"))
tissue.motifs <- tissue.motifs[which( ! tissue.motifs %in% c("SRF"))]

# remove tissue motifs in rhyth
rhyth.motifs <- rhyth.motifs[which(!rhyth.motifs %in% intersect(rhyth.motifs, tissue.motifs))]

# cross prods
mat.rhyth3 <- subset(mat.fgbg.3, select = intersect(rhyth.motifs, colnames(mat.fgbg.3)))
mat.tiss3 <- subset(mat.fgbg.3, select = intersect(tissue.motifs, colnames(mat.fgbg.3)))
mat.rhythtiss3 <- CrossProductTwoSets(mat.rhyth3, mat.tiss3)

mat.fgbg.cross.rhythtiss3 <- cbind(mat.fgbg.3, mat.rhythtiss3)
# remove columns with 0 variance
mat.fgbg.cross.rhythtiss3[which(colSums(mat.fgbg.cross.rhythtiss3) == 0)] <- list(NULL)

jlambda <- 0.035  # liv only
out.cross.rhythtiss3 <- PenalizedLDA(mat.fgbg.cross.rhythtiss3, labels3, lambda = jlambda, K = 2, standardized = FALSE)

# plot pretty
vec.length <- sqrt(out.cross.rhythtiss3$discrim[, 1]^2 + out.cross.rhythtiss3$discrim[, 2]^2)

jsize.cutoff <- 0.1
jsize.pairs.cut <- sapply(vec.length, function(jsize){
  if (jsize > jsize.cutoff){
    return(jsize)
  } else {
    return(0)
  }
})

labels <- names(out.cross.rhythtiss3$x)
labels.cut <- mapply(function(jlab, jsize){
  if (jsize <= 0){
    return("")
  } else {
    return(jlab)
  }
}, labels, jsize.pairs.cut)

dat.plot <- data.frame(x = out.cross.rhythtiss3$discrim[, 1],
                       y = out.cross.rhythtiss3$discrim[, 2],
                       motif = labels.cut,
                       vec.length = vec.length,
                       vec.length.cut = jsize.pairs.cut)
dat.labs <- subset(dat.plot, vec.length.cut > 0)

m <- ggplot(dat.plot, aes(x = x, y = y)) + 
  geom_point(size = 0.01) + 
  geom_text_repel(data = dat.labs, aes(x = x, y = y, label = motif), size = 2.5) + 
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + theme_bw() + 
  theme(aspect.ratio = 0.33, legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  xlab("Motif loadings separating liver and kidney DHS peaks") + ylab("Motif loadings separating rhythmic and flat DHS peaks")
print(m)


# Remove redundancies -----------------------------------------------------

fits.liv <- subset(fits.bytiss, tissue == "Liver_SV129")
amp.hash <- hash(as.character(fits.liv$gene), fits.liv$amp)
phase.hash <- hash(as.character(fits.liv$gene), fits.liv$phase)
pval.hash <- hash(as.character(fits.liv$gene), fits.liv$pval)

discrim <- out.cross.rhythtiss3$discrim
rownames(discrim) <- names(out.cross.rhythtiss3$x)
# discrim.hits <- discrim  # take all
# discrim.hits <- discrim[which(apply(discrim, 1, min) > 0.0125), ]  # must separate both liver and rhythms
# discrim.hits <- discrim[which(apply(discrim, 1, function(row) ifelse (row[[1]] > 0.01 & row[[2]] > 0.1, TRUE, FALSE)) == TRUE), ]  # must separate both liver and rhythms
# discrim.hits <- discrim[grepl("RORA;|SRF;|FOX\.F1;", rownames(discrim)), ]
# discrim.hits <- discrim[grepl("RORA;|SRF;", rownames(discrim)), ]
# discrim.hits <- discrim[grepl("RORA;|SRF;|FOX\\.F1\\.F2\\.J1\\.;", rownames(discrim)), ]
# discrim.hits <- discrim[grepl("RORA;|SRF;|bHLH_family;", rownames(discrim)), ]
# discrim.hits <- discrim[grepl("RORA;|SRF;", rownames(discrim)), ]
discrim.hits <- discrim[grepl("RORA;|bHLH_family;|SRF;", rownames(discrim)), ]
# discrim.hits <- discrim[grepl("FOX\\.F1\\.F2\\.J1\\.;", rownames(discrim)), ]
motifs.hits <- rownames(discrim.hits)
print("Motif hits: ")
print(motifs.hits)
if (length(motifs.hits) == 0) stop("No motifs!")
# add RORA and SRF
# motifs.hits <- c(motif.hits, c("RORA", "SRF"))

mat.M <- mat.fgbg.cross.rhythtiss3[which(labels3 == 1), motifs.hits]
# Plot loadings
# FILTER MOTIFS OPTIONALLY
# mat.M <- mat.M[, grepl(";", colnames(mat.M))]  # interactions only
# mat.M <- mat.M[, grepl(";", colnames(mat.M))]  # interactions only
# mat.M <- mat.M[, !grepl("FOX\\.F1\\.F2", colnames(mat.M))]  # no FOX
mat.M <- t(scale(t(mat.M), center = FALSE, scale = FALSE))
mat.M <- mat.M[which(rowSums(mat.M) > 0), ]
mat.M <- mat.fgbg.cross.rhythtiss3[which(labels3 == 1), motifs.hits]

# Do for penalized
library(PMA)
# K <- length(motifs.hits)

mat.M <- as.matrix(mat.M)
mat.pmd <- PMD(mat.M, type = c("standard"), sumabs = 0.4, center=FALSE, rnames = rownames(mat.M), cnames = colnames(mat.M), K=K)
print(mat.pmd)

rownames(mat.pmd$u) <- mat.pmd$rnames
rownames(mat.pmd$v) <- mat.pmd$cnames

pc1 <- 1
pc2 <- 3
plot(x = mat.pmd$v[, pc1], y = mat.pmd$v[, pc2])
text(x = mat.pmd$v[, pc1], y = mat.pmd$v[, pc2], labels = mat.pmd$cnames)

# map onto the two eigenpeak space
# mat.M.orth <- mat.M %*% mat.pmd$v  %*% diag(mat.pmd$d)  # since v's are all negative, flip to negative 1 it should be equivalent. 
# mat.M.orth <- mat.M %*% mat.pmd$v
mat.pmd$v <- sweep(mat.pmd$v, MARGIN = 2, STATS = sign(colSums(mat.pmd$v)), FUN = "*")  # fix sign so mean is positive
mat.M.orth <- mat.M %*% mat.pmd$v  # just rotate by orthonormal vectors
# mat.M.orth <- mat.pmd$u %*% diag(mat.pmd$d) %*% t(mat.pmd$v)
# mat.M.orth <- abs(mat.M.orth)  # make positive
colnames(mat.M.orth) <- paste("motifs_pc", seq(ncol(mat.pmd$v)), sep = "_")

N.peaks <- melt(mat.M.orth, varnames = c("peakgene", "motif"), value.name = "sitecount")
N.peaks$gene <- as.factor(sapply(as.character(N.peaks$peakgene), function(p) strsplit(p, ";")[[1]][[2]]))

detach("package:PMA", unload=TRUE)
detach("package:plyr", unload=TRUE)
library(dplyr)
# or else dplyr has bugs due to plyr
N.gene <- N.peaks %>%
  group_by(gene, motif) %>%
  summarise(sitecount = sum(sitecount))

# remove outliers!?!?!?
# outliers <- c("Pde9a", "Slc4a4", "Slc38a4", "Lhpp", "Upp2", "Kif1b")
# outliers <- c("Pde9a")
# N.gene <- subset(N.gene, !gene %in% outliers)

# from write_N_sitecounts_table_for_mara.R
mat.liver.cross <- dcast(N.gene, formula = gene ~ motif, value.var = "sitecount", fun.aggregate = sum, fill = 0)
rownames(mat.liver.cross) <- mat.liver.cross$gene; mat.liver.cross$gene <- NULL
mat.liver.cross <- mat.liver.cross[which(apply(mat.liver.cross, 1, sum) != 0), ]
mat.liver.cross <- cbind(data.frame(gene = rownames(mat.liver.cross)), mat.liver.cross)

outdir <- "/home/yeung/projects/tissue-specificity/data/sitecounts/motevo/liver_kidney_sitecounts_tissuespecpeaks_cutofflow0_pma_decorrelate"
dir.create(outdir)

outf.mat.cross <- file.path(outdir, paste0("2D.posterior.multigene.distfilt.morenonliv.bugfixed.liverWTKO.40000.cutoff.3.cutofflow0.method.g1001.", descrip, ".K_", K, ".mat"))
print("Writing to:")
print(outf.mat.cross)
write.table(mat.liver.cross, file = outf.mat.cross,
            quote = FALSE, sep = "\t", row.names = FALSE)
# 
# # downstream analysis of phases of these 141 genes
# mat.liver.cross.ds <- mat.liver.cross
# mat.liver.cross.ds$phase <- sapply(as.character(mat.liver.cross$gene), function(g) phase.atg[[g]])
# mat.liver.cross.ds$amp <- sapply(as.character(mat.liver.cross$gene), function(g) amp.atg[[g]])
# 
# pc1 <- "motifs_pc_1"
# pc2 <- "motifs_pc_2"
# ggplot(mat.liver.cross.ds, aes_string(x = "phase", y = pc1)) + geom_point() + ggtitle("PC1 motifs")
# ggplot(mat.liver.cross.ds, aes_string(x = "phase", y = pc2)) + geom_point() + ggtitle("PC2 motifs")
# 
# RUN MARA
# cd /home/yeung/projects/tissue-specificity/scripts/liver_kidney_WTKO
# TRY multiple mat (deprecated)
# bash /home/yeung/projects/tissue-specificity/scripts/liver_kidney_WTKO/run_run_filter_mara_pipeline.promoters.bygenelist.sh /home/yeung/projects/tissue-specificity/data/sitecounts/motevo/liver_kidney_sitecounts_tissuespecpeaks_cutofflow0_cross_only_pma_decorrelate /home/yeung/projects/tissue-specificity/results/MARA.liver_kidney/cross_only_decorrelated
# Single MAT 
# bash /home/yeung/projects/tissue-specificity/scripts/liver_kidney_WTKO/run_run_filter_mara_pipeline.promoters.bygenelist.singlemat.sh /home/yeung/projects/tissue-specificity/data/sitecounts/motevo/liver_kidney_sitecounts_tissuespecpeaks_cutofflow0_cross_only_pma_decorrelate /home/yeung/projects/tissue-specificity/results/MARA.liver_kidney/cross_only_decorrelated

marascript <- "/home/yeung/projects/tissue-specificity/scripts/liver_kidney_WTKO/run_run_filter_mara_pipeline.promoters.bygenelist.singlemat.sh"
nmat <- outf.mat.cross
outmain <- paste0("/home/yeung/projects/tissue-specificity/results/MARA.liver_kidney/mara_out.", descrip, ".K_", K)
cmd <- paste("bash", marascript, nmat, outmain)
system(cmd)

# LOAD RESULTS
# indir <- "/home/yeung/projects/tissue-specificity/results/MARA.liver_kidney/cross_only_decorrelated/2D.posterior.multigene.distfilt.morenonliv.bugfixed.liverWTKO.40000.cutoff.3.cutofflow0.method.g1001.mat/atger_with_kidney"
# indir <- "/home/yeung/projects/tissue-specificity/results/MARA.liver_kidney/cross_only_decorrelated/2D.posterior.multigene.distfilt.morenonliv.bugfixed.liverWTKO.40000.cutoff.3.cutofflow0.method.g1001.mat/atger_with_kidney.bugfixed
# indir <- "/home/yeung/projects/tissue-specificity/results/MARA.liver_kidney/cross_only_decorrelated.bugfixed/2D.posterior.multigene.distfilt.morenonliv.bugfixed.liverWTKO.40000.cutoff.3.cutofflow0.method.g1001.mat/atger_with_kidney.bugfixed"
# indir <- "/home/yeung/projects/tissue-specificity/results/MARA.liver_kidney/cross_only_decorrelated.Pde9aremoved/2D.posterior.multigene.distfilt.morenonliv.bugfixed.liverWTKO.40000.cutoff.3.cutofflow0.method.g1001.Pde9a.removed.mat/atger_with_kidney.bugfixed"
# indir <- "/home/yeung/projects/tissue-specificity/results/MARA.liver_kidney/cross_only_decorrelated.Pde9.Slc4a4.removed/2D.posterior.multigene.distfilt.morenonliv.bugfixed.liverWTKO.40000.cutoff.3.cutofflow0.method.g1001.Pde9a.Slc4a4.removed.mat/atger_with_kidney.bugfixed"
# indir <- "/home/yeung/projects/tissue-specificity/results/MARA.liver_kidney/cross_only_decorrelated.Pde9.Slc4a4.removed/2D.posterior.multigene.distfilt.morenonliv.bugfixed.liverWTKO.40000.cutoff.3.cutofflow0.method.g1001.Pde9a.Slc4a4.Slc38a4.removed.mat/atger_with_kidney.bugfixed"
# indir <- "/home/yeung/projects/tissue-specificity/results/MARA.liver_kidney/cross_only_decorrelated.all_motifs_decorrelate.removed/2D.posterior.multigene.distfilt.morenonliv.bugfixed.liverWTKO.40000.cutoff.3.cutofflow0.method.g1001.Pde9a.Slc4a4.Slc38a4.removed.mat/atger_with_kidney.bugfixed"
# indir <- "/home/yeung/projects/tissue-specificity/results/MARA.liver_kidney/cross_only_decorrelated.all_motifs_decorrelate.removed/2D.posterior.multigene.distfilt.morenonliv.bugfixed.liverWTKO.40000.cutoff.3.cutofflow0.method.g1001.allmotifsdecorrlate.removed.mat/atger_with_kidney.bugfixed"
# indir <- "/home/yeung/projects/tissue-specificity/results/MARA.liver_kidney/cross_only_decorrelated.all_motifs_decorrelate.removed.K4/2D.posterior.multigene.distfilt.morenonliv.bugfixed.liverWTKO.40000.cutoff.3.cutofflow0.method.g1001.K4.mat/atger_with_kidney.bugfixed"
# indir <- "/home/yeung/projects/tissue-specificity/results/MARA.liver_kidney/cross_only_decorrelated.K_10/atger_with_kidney.bugfixed"
# indir <- "/home/yeung/projects/tissue-specificity/results/MARA.liver_kidney/many_motifs_decorrelated.nofiltergenes/atger_with_kidney.bugfixed"
# indir <- "/home/yeung/projects/tissue-specificity/results/MARA.liver_kidney/many_motifs_decorrelated.nofiltergenes.noabs/atger_with_kidney.bugfixed"
inmain <- outmain
indir <- file.path(inmain, "atger_with_kidney.bugfixed")
# indir <- "/home/yeung/projects/tissue-specificity/results/MARA.liver_kidney/many_motifs_decorrelated.nofiltergenes.noabs/atger_with_kidney.bugfixed"
act.long <- LoadActivitiesLongKidneyLiver(indir, collapse.geno.tissue = TRUE, shorten.motif.name = FALSE)
# rename motifs based on the motifs with non-zero entries
act.long$pc <- sapply(as.character(act.long$gene), function(g) as.numeric(strsplit(g, "_")[[1]][[3]]), USE.NAMES = FALSE)
motif.decor <- sapply(act.long$pc, function(p) PcToMotif(mat.pmd$v, p), USE.NAMES = FALSE)
act.long$gene <- motif.decor


omega <- 2 * pi / 24
act.complex <- act.long %>%
  group_by(gene, tissue) %>%
  do(ProjectToFrequency2(., omega, add.tissue=TRUE))


s.act <- SvdOnComplex(act.complex, value.var = "exprs.transformed")

jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)
# jtitle <- gsub(pattern = "\\.", replacement = "\n", basename(indirmain))

max.labs <- 2
jtitle <- ""
for (comp in seq(1)){
  eigens.act <- GetEigens(s.act, period = 24, comp = comp, adj.mag = TRUE, constant.amp = 4, label.n = K, jtitle = jtitle)
  print(eigens.act$u.plot)
}

# # check what the motifs are
pc1 <- 4
pc2 <- 3

library(wordcloud)
wordcloud::textplot(x = mat.pmd$v[, pc1], y = mat.pmd$v[, pc2], words = mat.pmd$cnames, xlab = paste("PC", pc1), ylab = paste("PC", pc2))
# pc2 <- 3
# # pc2 <- 2
# # pc2 <- 4
# # pc1 <- 10
# # pc1 <- 7
# 
# ggplot(subset(act.long, gene == "motifs_pc_5"), aes(x = time, y = exprs)) + geom_point() + geom_line() + facet_wrap(~tissue)
# 
# ggplot(subset(act.long, gene == "motifs_pc_2"), aes(x = time, y = exprs)) + geom_point() + geom_line() + facet_wrap(~tissue)
# ggplot(subset(act.long, gene == "motifs_pc_3"), aes(x = time, y = exprs)) + geom_point() + geom_line() + facet_wrap(~tissue)
# 
# 
# # Compare with gene expression --------------------------------------------
# 
# load("Robjs/liver_kidney_atger_nestle/dat.freq.bugfixed.Robj")
# i <- 1
# 
# genes.filt <- as.character(mat.liver.cross$gene)
# s <- SvdOnComplex(subset(dat.freq, gene %in% genes.filt), value.var = "exprs.transformed")
# eigens <- GetEigens(s, period = 24, comp = i, label.n = 30, eigenval = TRUE, adj.mag = TRUE, constant.amp = 4, peak.to.trough = TRUE)
# jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)
# multiplot(eigens$u.plot, eigens$v.plot, layout = jlayout)
