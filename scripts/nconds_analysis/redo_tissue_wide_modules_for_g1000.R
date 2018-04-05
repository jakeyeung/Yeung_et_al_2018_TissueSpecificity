# 2016-07-27
# Jake Yeung
# Since Liver and Kidney was done on g=1000, redo eerything for g=1000

rm(list=ls())

source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/NcondsFunctions.R")
source("scripts/functions/SvdFunctions.R")

load("/home/yeung/projects/tissue-specificity/Robjs/nconds_g1000_11_tissues/fits_long.11_tiss_3_max.g1000.bestmodel.Robj", v=T)
load("Robjs/dat.long.fixed_rik_genes.Robj", v=T)
load("Robjs/fits.relamp.Robj", v=T)
load("Robjs/dat.complex.fixed_rik_genes.Robj", v=T)
load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.15.phase_sd_maxdiff_avg.Robj", v=T)

filter.by.amp <- 0.15

# annotate ----------------------------------------------------------------

if (filter.by.amp){
  fits.long$model <- mapply(FilterModelByAmp, fits.long$model, fits.long$param.list, MoreArgs = list(amp.cutoff = 0.15))
}
fits.long$n.params <- sapply(fits.long$model, function(m) return(length(strsplit(as.character(m), ";")[[1]])))
fits.long$n.rhyth <- sapply(fits.long$model, GetNrhythFromModel)
fits.long$amp.avg <- mapply(GetAvgAmpFromParams, fits.long$param.list, fits.long$model)
fits.long$phase.sd <- mapply(GetSdPhaseFromParams, fits.long$param.list, fits.long$model)
fits.long$phase.maxdiff <- mapply(GetMaxPhaseDiffFromParams, fits.long$param.list, fits.long$model)
fits.long$phase.avg <- mapply(GetPhaseFromParams, fits.long$param.list, fits.long$model)

if (filter.by.amp == 0.15){
  save(fits.long, file = "Robjs/nconds_g1000_11_tissues/fits_long.11_tiss_3_max.g1000.bestmodel.filteramp.0.15.Robj")
}

# Sanity check ------------------------------------------------------------

genes.tw <- as.character(subset(fits.best, n.rhyth >= 8)$gene)  # to compare
genes.tw <- as.character(subset(fits.long, n.rhyth >= 8)$gene)

s.tw <- SvdOnComplex(subset(dat.complex, gene %in% genes.tw), value.var = "exprs.transformed")
eigens.tw <- GetEigens(s.tw, period = 24, comp = 1, label.n = 15, eigenval = TRUE, adj.mag = TRUE, constant.amp = 4)
jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)
multiplot(eigens.tw$u.plot, eigens.tw$v.plot, layout = jlayout)


# Write gene expression directory -----------------------------------------



# Run MARA  ---------------------------------------------------------------

Npath <- "/home/yeung/projects/tissue-specificity/data/sitecounts/motevo/sitecount_matrix_geneids"
Npath.outmain <- "/home/yeung/projects/tissue-specificity/data/sitecounts/motevo/promoters_hogenesch"
dir.create(Npath.outmain)
Npath.out <- paste0(Npath.outmain, "sitecount.tissuewide.filteramp.", filter.by.amp, ".mat")
N <- read.table(Npath, header=TRUE)

N.sub <- subset(N, Gene.ID %in% genes.tw)
print(paste("Writing", nrow(N.sub), "promoters to file:", Npath.out))

write.table(N.sub, file = Npath.out, append = FALSE, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)

# Run MARA ----------------------------------------------------------------

marascript <- "/home/yeung/projects/tissue-specificity/scripts/liver_kidney_WTKO/run_run_filter_mara_pipeline.promoters.bygenelist.singlemat.sh"
geneexprsdir <- "/home/yeung/projects/tissue-specificity/data/gene_lists/expressed_genes_deseq_int/expressed_genes_deseq_int.centeredTRUE"
nmat <- Npath.out

outbase <- "/home/yeung/projects/tissue-specificity/results/MARA.hogenesch"
dir.create(outbase)
outmain <- file.path(outbase, paste0("promoters.tissuewide.filteramp.", filter.by.amp, ".mat"))
# outmain <- paste0("/home/yeung/projects/tissue-specificity/results/MARA.liver_kidney/promoters.", jmod)

cmd <- paste("bash", marascript, nmat, outmain, geneexprsdir)
if (!dir.exists(outmain)){
  system(cmd)
}

# Downstream --------------------------------------------------------------

indir <- file.path(outmain, "expressed_genes_deseq_int.centeredTRUE")
source("scripts/functions/LoadActivitiesLong.R")
act.long <- LoadActivitiesLong(indir)
# rename motifs based on the motifs with non-zero entries

omega <- 2 * pi / 24
act.complex <- subset(act.long, tissue != "WFAT") %>%
  group_by(gene, tissue) %>%
  do(ProjectToFrequency2(., omega, add.tissue=TRUE))

s.act <- SvdOnComplex(act.complex, value.var = "exprs.transformed")


jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)
# jtitle <- gsub(pattern = "\\.", replacement = "\n", basename(indirmain))

max.labs <- 20
jtitle <- ""
for (comp in seq(1)){
  eigens.act <- GetEigens(s.act, period = 24, comp = comp, adj.mag = TRUE, constant.amp = 4, label.n = max.labs, jtitle = jtitle, peak.to.trough = TRUE)
  print(eigens.act$u.plot)
  multiplot(eigens.act$u.plot, eigens.act$v.plot, cols = 2)
}
