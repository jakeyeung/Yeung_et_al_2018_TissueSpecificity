# 2016-10-26
# Jake YEung
# go_terms_analysis_with_phase_amplitude.R
# copied from liver kidney WT KO analysis
# Create a phase-specific annotation (sliding window?) of functional enrichment. 

rm(list=ls())
start <- Sys.time()

textsize <- 7
dotsize <- 2
remove.kidney.outliers <- TRUE

library(ggplot2)
# detach("package:PMA", unload=TRUE)
# detach("package:plyr", unload=TRUE)
# detach("package:ggplot2", unload=TRUE)
# detach("package:reshape2", unload=TRUE)
library(dplyr)
library(parallel)

setwd("/home/yeung/projects/tissue-specificity")

source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/NcondsFunctions.R")
source("scripts/functions/SvdFunctions.R")
source("scripts/functions/LoadActivitiesLong.R")
source("scripts/functions/LiverKidneyFunctions.R")
source("scripts/functions/PlotActivitiesFunctions.R")
source("scripts/functions/FourierFunctions.R")
source("scripts/functions/GetTFs.R")
source("scripts/functions/LdaFunctions.R")
source("scripts/functions/HandleMotifNames.R")
source("scripts/functions/RemoveP2Name.R")
source("scripts/functions/GetTopMotifs.R")
source("scripts/functions/NcondsAnalysisFunctions.R")
source("scripts/functions/ModelStrToModel.R")
source("scripts/functions/ListFunctions.R")
source("scripts/functions/AnalyzeGeneEnrichment.R")
source("scripts/functions/BiomartFunctions.R")


WriteListToFile <- function(lst, outf){
  sink(file = outf)
  for (g in lst){
    cat(g)
    cat("\n")
  }
  sink()
  return(NA)
}


Vectorize(IsBtwnTimes <- function(phase, tstart, tend){
  # check if phase (between 0 to 24, is between tstart and tend, considering the modulo)
  if (tend > tstart){
    # easy case
    is.btwn <- phase >= tstart & phase <= tend
  } else {
    # harder case, must consider the modulo
    is.btwn <- phase >= tstart | phase <= tend
  }
  # replace NAs with FALSE
  is.btwn[which(is.na(is.btwn))] <- FALSE
  return(is.btwn)
}, vectorize.args="phase")

load("Robjs/dat.long.fixed_rik_genes.Robj", v=T)  # hogenesch

# Annotate each model with its major GO term ------------------------------

jmeth <- "g=1001"
load("Robjs/liver_kidney_atger_nestle/dat.freq.bugfixed.Robj", v=T)  # LivKid
load("Robjs/liver_kidney_atger_nestle/fits.long.multimethod.filtbest.staggeredtimepts.bugfixed.annotated.Robj", v=T)
fits.long.filt <- subset(fits.long.filt, method == jmeth)
load("Robjs/liver_kidney_atger_nestle/dat.long.liverkidneyWTKO.bugfixed.Robj", v=T); dat.wtko <- dat.long; rm(dat.long)
load("Robjs/liver_kidney_atger_nestle/fits.bytiss.bugfixed.Robj", v=T)

fits.long.filt$n.params <- sapply(fits.long.filt$model, function(m) return(length(strsplit(as.character(m), ";")[[1]])))
fits.long.filt$n.rhyth <- sapply(fits.long.filt$model, GetNrhythFromModel)

# Remove kidney outliers (optional)
if (remove.kidney.outliers){
  # Kidney_SV129 genes contain some weird outliers, remove them
  outliers <- c("Cox8b", "Ucp1", "Cidea", "Flg", "Nr4a2")
  fits.long.filt <- subset(fits.long.filt, !gene %in% outliers)
  dat.wtko <- subset(dat.wtko, !gene %in% outliers)
}
dat.wtko <- StaggeredTimepointsLivKid(dat.wtko)
dat.wtko.collapsed <- CollapseTissueGeno(dat.wtko)

# Plot for all genes in module --------------------------------------------

# from tissue_specificity_paper3.R
jmod1 <- "Liver_SV129,Kidney_SV129,Liver_BmalKO,Kidney_BmalKO-Liver_SV129,Liver_BmalKO.Kidney_SV129,Kidney_BmalKO"
jmod2 <- "Liver_SV129,Kidney_SV129.Liver_BmalKO,Kidney_BmalKO-Liver_SV129,Kidney_SV129"
jtiss.lst <- list(c(ModelStrToModel(jmod1), "MF"),
                  c(ModelStrToModel(jmod2), "BP"),
                  c("Liver_SV129", "BP"), 
                  c("Liver_SV129,Liver_BmalKO", "BP"), 
                  c("Kidney_SV129", "MF"), 
                  c("Kidney_SV129,Kidney_BmalKO", "MF"), 
                  c("Liver_BmalKO", "BP"))

plots <- expandingList()
jmod.long <- "Liver_SV129,Liver_BmalKO"
jonto <- "BP"
comp <- 1
jcutoff <- 0.5
jtiss.onto <- paste(jmod.long, jonto)

genes <- as.character(subset(fits.long.filt, model %in% jmod.long)$gene)
dat.sub <- subset(dat.freq, gene %in% genes)
s <- SvdOnComplex(dat.sub, value.var = "exprs.transformed")
eigens <- GetEigens(s, period = 24, comp = comp, label.n = 25, eigenval = TRUE, adj.mag = TRUE, constant.amp = textsize, peak.to.trough = TRUE, label.gene = c("Egr1", "Mafb", "Tfcp2"))

plots$add(eigens$u.plot + ylab("ZT") + ggtitle(""))
plots$add(eigens$v.plot + ylab("ZT") + xlab("Tissue Weights") + ggtitle(""))

# genes.bg <- as.character(subset(fits.long.filt, model == "")$gene)
genes.bg <- as.character(fits.long.filt$gene)
# genes.bg <- as.character(subset(fits.long.filt, model %in% jmod.long & phase.avg > 12 & phase.avg < 24)$gene)
# genes.bg.ensembl <- Gene2Ensembl(genes.bg)
# genes.fg <- as.character(subset(fits.long.filt, model %in% jmod.long & phase.avg > 0 & phase.avg < 12)$gene)
# genes.fg.ensembl <- Gene2Ensembl(genes.fg)

# write fg and bg to separate file to load into David:
genedir <- "/home/yeung/projects/tissue-specificity/data/gene_lists/GO_analysis"

# GO terms for Liver WTKO module
GOterms <- c("GO:0006260", "GO:0042254", "GO:0032868", "GO:0043434", "GO:0046326")
# plot in 6 hour intervals
for (tstart in c(3, 9, 21)){
  tend <- tstart + 6
  if (tend > 24){
    tend <- tend - 24
  }
  print(paste(tstart, tend))
  genes <- as.character(subset(fits.long.filt, model %in% jmod.long & IsBtwnTimes(phase.avg, tstart, tend))$gene)
  print(paste("Ngenes", length(genes)))
  genes.ensembl <- Gene2Ensembl(genes)
  enrichment <- GetGOEnrichment(genes.bg, genes, fdr.cutoff = 1, ontology = jonto, show.top.n = Inf)
  # outf <- file.path(genedir, paste0(gsub(",", "-", jmod.long), ".", tstart, "to", tend, ".txt"))
  # WriteListToFile(genes.ensembl, outf)
}
# 
# sink(file = file.path(genedir, paste0(gsub(",", "-", jmod.long), ".12to24.bg.txt")))
# for (g in genes.bg.ensembl){
#   cat(g)
#   cat("\n")
# }
# sink()
# 
# enrichment <- GetGOEnrichment(genes.bg, genes.fg, fdr.cutoff = 1, ontology = jonto, show.top.n = 8)
# 
# m.go <- ggplot(enrichment, aes(x = Term, y = minuslogpval)) + geom_bar(stat = "identity") +
#   ylab("-log10(P-value), Fisher's exact test") +
#   xlab("") +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1), aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   ggtitle(jtiss.onto)
# m.go
# 
# m <- PlotGeneModuleWithGO(dat.sub, enrichment, jtitle = paste(jtiss.onto, collapse = "\n"), text.size = textsize, dot.size = dotsize, comp = comp, legend.pos = "none", label.gene = c("Egr1", "Mafb", "Tfcp2"), label.GO.terms.only = TRUE, unlabeled.alpha = 0.3)
# m2 <- PlotGeneModuleWithGO(dat.sub, enrichment, jtitle = paste(jtiss.onto, collapse = "\n"), text.size = textsize, dot.size = dotsize, comp = comp, legend.pos = "right", label.gene = c("Egr1", "Mafb", "Tfcp2"), label.GO.terms.only = TRUE, unlabeled.alpha = 0.3)
# plots$add(m.go)
# plots$add(m + ylab("ZT") + ggtitle(""))
# plots$add(m2 + ylab("ZT") + ggtitle(""))
# 
# plotdir <- "/home/yeung/projects/tissue-specificity/plots/liver_kidney_modules_with_GO"
# pdf(file.path(plotdir, paste0("liv_kid_with_GO.rm_outliers.", remove.kidney.outliers, ".pdf")))
#   mclapply(jtiss.lst, function(jtiss.onto){
#     # jtiss.onto <- jtiss.lst[[1]]  # c(jmodels, ontology), remove last element to get jmodels, last element is ontology
#     source("scripts/functions/AnalyzeGeneEnrichment.R")
#     jcutoff <- 0.5  # fdr cutoff for enrichment
#     comp <- 1
#     plots <- expandingList()
#     jmod.long <- jtiss.onto[-length(jtiss.onto)]
#     jonto <- jtiss.onto[length(jtiss.onto)]
#     
#     genes <- as.character(subset(fits.long.filt, model %in% jmod.long)$gene)
#     dat.sub <- subset(dat.freq, gene %in% genes)
#     s <- SvdOnComplex(dat.sub, value.var = "exprs.transformed")
#     eigens <- GetEigens(s, period = 24, comp = comp, label.n = 25, eigenval = TRUE, adj.mag = TRUE, constant.amp = textsize, peak.to.trough = TRUE, label.gene = c("Egr1", "Mafb", "Tfcp2"))
#     
#     plots$add(eigens$u.plot + ylab("ZT") + ggtitle(""))
#     plots$add(eigens$v.plot + ylab("ZT") + xlab("Tissue Weights") + ggtitle(""))
#     
#     genes.bg <- as.character(subset(fits.long.filt)$gene)
#     genes.fg <- as.character(subset(fits.long.filt, model %in% jmod.long)$gene)
#     enrichment <- GetGOEnrichment(genes.bg, genes.fg, fdr.cutoff = jcutoff, ontology = jonto, show.top.n = 8)
#     # merge GO terms (manual for the moment)
#     if (jtiss.onto[[1]] == "Liver_SV129,Kidney_SV129;Liver_BmalKO,Kidney_BmalKO"){
#       # merge circ reg and rhyth process into one term
#       go.terms <- c("GO:0032922", "GO:0048511")
#       enrichment <- MergeGOTerms(enrichment, go.terms, new.go.term = "circadian or rhythmic process")
#     }
#     # plot enrichment GO terms
#     m.go <- ggplot(enrichment, aes(x = Term, y = minuslogpval)) + geom_bar(stat = "identity") +
#       ylab("-log10(P-value), Fisher's exact test") +
#       xlab("") +
#       theme_bw() +
#       theme(axis.text.x = element_text(angle = 45, hjust = 1), aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#       ggtitle(jtiss.onto)
#     m <- PlotGeneModuleWithGO(dat.sub, enrichment, jtitle = paste(jtiss.onto, collapse = "\n"), text.size = textsize, dot.size = dotsize, comp = comp, legend.pos = "none", label.gene = c("Egr1", "Mafb", "Tfcp2"), label.GO.terms.only = TRUE, unlabeled.alpha = 0.3)
#     m2 <- PlotGeneModuleWithGO(dat.sub, enrichment, jtitle = paste(jtiss.onto, collapse = "\n"), text.size = textsize, dot.size = dotsize, comp = comp, legend.pos = "right", label.gene = c("Egr1", "Mafb", "Tfcp2"), label.GO.terms.only = TRUE, unlabeled.alpha = 0.3)
#     plots$add(m.go)
#     plots$add(m + ylab("ZT") + ggtitle(""))
#     plots$add(m2 + ylab("ZT") + ggtitle(""))
#     return(plots$as.list())
#   }, mc.cores = length(jtiss.lst))
# dev.off()


# Check Kidney_SV129 genes for false positives ----------------------------

# jmod <- "Kidney_SV129,Kidney_BmalKO"
# fits.sub <- subset(fits.long.filt, model == jmod)
# fits.sub <- fits.sub[order(fits.sub$amp.avg, decreasing = TRUE), ]
# # order from highest amp to lowest amp
# jgenes <- as.character(fits.sub$gene)
# 
# pdf(file.path(plotdir, paste0(jmod, "_sanity_check.pdf")))
# for (g in jgenes){
#   print(PlotGeneTissuesWTKO(subset(dat.wtko, gene == g), jtitle = g))
# }
# dev.off()

# 
# 
# # Highlight genes assigned to GO terms ------------------------------------
# genes.bg <- as.character(subset(fits.long.filt)$gene)
# genes.fg <- as.character(subset(fits.long.filt, model %in% jmod.long)$gene)
# enrichment <- GetGOEnrichment(genes.bg, genes.fg, fdr.cutoff = jcutoff, ontology = jonto, show.top.n = 10)
# 
# m1 <- ggplot(subset(enrichment, !is.na(Term)), aes(x = Term, y = minuslogpval)) + geom_bar(stat = "identity") + 
#   ylab("-log10(P-value), Fisher's exact test") + 
#   xlab("") +
#   theme_bw() + 
#   theme(axis.text.x = element_text(angle = 45, hjust = 1), aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   ggtitle(jmod.long)
# print(m1)
# 
# # annotate dat.freq with GO terms for each gene (take most significant enrichment)
# 
# # label clocks to based on GOterm
# go.hash <- hash()
# for (i in seq(nrow(enrichment))){
#   genes.vec <- enrichment$genes[[i]]
#   term <- as.character(enrichment$Term[[i]])
#   for (g in genes.vec){
#     if (is.null(go.hash[[g]])){
#       go.hash[[g]] <- term
#     } 
#   }
# }
# 
# dat.sub <- subset(dat.freq, gene %in% genes.fg)
# s.sub <- SvdOnComplex(dat.sub, value.var = "exprs.transformed")
# 
# # label dat.freq with goterm
# dat.sub$term <- sapply(as.character(dat.sub$gene), function(g){
#   if (!is.null(go.hash[[g]])){
#     return(go.hash[[g]])
#   } else {
#     return(NA)
#   }
# })
# 
# eig <- GetEigens(s.sub, period = 24, comp = comp, label.n = 25, eigenval = TRUE, adj.mag = TRUE, constant.amp = 4, peak.to.trough = TRUE)
# 
# omega <- 2 * pi / 24
# ampscale <- 2
# vec.complex <- eig$eigensamp 
# labels <- names(vec.complex)
# 
# dat <- data.frame(amp = Mod(vec.complex) * ampscale,
#                  phase = ConvertArgToPhase(Arg(vec.complex), omega = omega),
#                  label = labels)
# 
# # add term BEFORE filtering out labels
# dat$term <- as.factor(sapply(as.character(dat$label), function(g){
#   if (g == ""){
#     return("")
#   }
#   if (!is.null(go.hash[[g]])){
#     return(go.hash[[g]])
#   } else {
#     return("")
#   }
# }))
# # make "" the first term for levels (so you get colours of dots wihtout the labels 
# if (any(dat$term == "")){
#   jterms <- as.character(enrichment$Term)
#   dat$term <- factor(as.character(dat$term), levels = c("", jterms))
# }
# 
# top.hits <- 25
# top.amps <- as.character(head(dat[order(dat$amp, decreasing = TRUE), ], n = top.hits)$label)
# dat$label <- sapply(as.character(dat$label), function(l) ifelse(l %in% top.amps, yes = l, no = ""))
# # label only top genes
# 
# 
# amp.max <- ceiling(max(dat$amp) * 2) / 2
# if (amp.max <= 1){
#   amp.step <- 0.5
# } else {
#   amp.step <- 1
# }
# 
# xlab <- ""
# ylab <- "Log2 Fold Change"
# constant.amp <- 6
# jtitle <- ""
# 
# cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# cbPalette <- cbPalette[1:length(levels(dat$term))]
# # change "" term to light gray
# cbPalette[levels(dat$term) == ""] <- "gray85"
# 
# m <- ggplot(data = dat, aes(x = amp, y = phase, label = label, colour = term)) + 
#   geom_point(size = 1) +
#   coord_polar(theta = "y") + 
#   xlab(xlab) +
#   ylab(ylab) +
#   ggtitle(jtitle) +
#   scale_y_continuous(limits = c(0, 24), breaks = seq(6, 24, 6)) + 
#   scale_x_continuous(limits = c(0, amp.max), breaks = seq(0, amp.max, length.out = 2)) + 
#   theme_bw() + 
#   geom_vline(xintercept = seq(0, amp.max, length.out = 2), colour = "grey50", size = 0.2, linetype = "dashed") +
#   geom_hline(yintercept = seq(6, 24, by = 6), colour = "grey50", size = 0.2, linetype = "solid") +
#   theme(panel.grid.major = element_line(size = 0.5, colour = "grey"), panel.grid.minor = element_blank(), 
#         panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position="bottom",
#         panel.border = element_blank(),
#         legend.key = element_blank(),
#         axis.ticks = element_blank(),
#         panel.grid  = element_blank())
# # expand_limits(x = 0)  # we want center to be 0
# 
# # add text
# df.txt <- subset(dat, label != "")
# if (constant.amp != FALSE){
#   # m <- m + geom_text_repel(data = df.txt, aes(x = amp, y = phase, label = label), size = constant.amp, colour = "black")
#   m <- m + geom_text_repel(data = df.txt, aes(x = amp, y = phase, label = label), size = constant.amp)
# } else {
#   # m <- m + geom_text_repel(data = df.txt, aes(x = amp, y = phase, size = amp, label = label), colour = "black")
#   m <- m + geom_text_repel(data = df.txt, aes(x = amp, y = phase, size = amp, label = label))
# }
# print(m + scale_colour_manual(values=cbPalette) + ggtitle(jonto))
# 
# # PlotEnrichmentGenes(dat.freq, enrichment, max.genes = 40, row.i = "max")
# 
