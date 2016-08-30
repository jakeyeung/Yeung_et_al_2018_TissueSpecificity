# 2016-08-30
# Jake YEung
# go_terms_analysis_with_phase_amplitude.R

rm(list=ls())
start <- Sys.time()

library(ggplot2)
# detach("package:PMA", unload=TRUE)
# detach("package:plyr", unload=TRUE)
# detach("package:ggplot2", unload=TRUE)
# detach("package:reshape2", unload=TRUE)
library(dplyr)

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

GetGOEnrichment <- function(genes.bg, genes.fg, fdr.cutoff, show.top.n = Inf, ontology="BP", wd = "/home/yeung/projects/tissue-specificity"){
  source(file.path(wd, "scripts/functions/AnalyzeGeneEnrichment.R"))
  enrichment <- AnalyzeGeneEnrichment(genes.bg, genes.fg, FDR.cutoff = 0.5, which.ontology = ontology, return.GOdata = TRUE)
  enrichment$minuslogpval <- -log10(as.numeric(enrichment$classicFisher))
  enrichment <- OrderDecreasing(enrichment, jfactor = "Term", jval = "minuslogpval")
  show.top.n.min <- min(nrow(enrichment), show.top.n)
  if (show.top.n.min == 0) return(NULL)
  enrichment <- enrichment[1:show.top.n.min, ]   # prevent taking more than you have enrichment
  # unload packages
  # sometimes topGO causes problems (unable to unload later), unload once you're done.
  detach(name = "package:topGO", unload = TRUE)
  detach(name = "package:org.Mm.eg.db", unload = TRUE)
  return(enrichment)
}

PlotEnrichmentGenes <- function(dat.freq, enrichment, max.genes, row.i = "max"){
  # take hit with most genes
  # if row.i is max, get row with most genes
  # otherwise use row.i as index for row
  if (row.i == "max"){
    row.i <- which(enrichment$Significant == max(enrichment$Significant, na.rm = TRUE))
  } else if (is.numeric(row.i)){
    row.i <- row.i
  } else {
    warning(paste("row.i is not max or numeric:", row.i))
  }
  go.genes <- enrichment$genes[[row.i]]
  go.term <- enrichment$Term[[row.i]]
  show.n.genes <- min(length(go.genes), max.genes)
  s <- SvdOnComplex(subset(dat.freq, gene %in% go.genes), value.var = "exprs.transformed")
  eigens <- GetEigens(s, period = 24, comp = comp, label.n = show.n.genes, eigenval = TRUE, adj.mag = TRUE, constant.amp = 4, peak.to.trough = TRUE, label.gene = FALSE)
  print(eigens$u.plot)
}

# Annotate each model with its major GO term ------------------------------

jmeth <- "g=1001"
load("Robjs/liver_kidney_atger_nestle/dat.freq.bugfixed.Robj", v=T)  # LivKid
load("Robjs/liver_kidney_atger_nestle/fits.long.multimethod.filtbest.staggeredtimepts.bugfixed.annotated.Robj", v=T)
fits.long.filt <- subset(fits.long.filt, method == jmeth)
load("Robjs/liver_kidney_atger_nestle/dat.long.liverkidneyWTKO.bugfixed.Robj", v=T); dat.wtko <- dat.long; rm(dat.long)
load("Robjs/liver_kidney_atger_nestle/fits.bytiss.bugfixed.Robj", v=T)

fits.long.filt$n.params <- sapply(fits.long.filt$model, function(m) return(length(strsplit(as.character(m), ";")[[1]])))
fits.long.filt$n.rhyth <- sapply(fits.long.filt$model, GetNrhythFromModel)

dat.wtko <- StaggeredTimepointsLivKid(dat.wtko)
dat.wtko.collapsed <- CollapseTissueGeno(dat.wtko)


# Plot for all genes in module --------------------------------------------

jmod1 <- "Liver_SV129,Kidney_SV129,Liver_BmalKO,Kidney_BmalKO-Liver_SV129,Liver_BmalKO.Kidney_SV129,Kidney_BmalKO"
jmod2 <- "Liver_SV129,Kidney_SV129.Liver_BmalKO,Kidney_BmalKO-Liver_SV129,Kidney_SV129"
jtiss.lst <- list(ModelStrToModel(jmod1),
                  ModelStrToModel(jmod2),
                  "Liver_SV129", 
                  "Liver_SV129,Liver_BmalKO", 
                  "Kidney_SV129", 
                  "Kidney_SV129,Kidney_BmalKO", 
                  "Liver_BmalKO")

jmod.long <- ModelStrToModel(jmod1)
jmod.long <- ModelStrToModel(jmod2)

jonto <- "BP"
jcutoff <- 0.5  # fdr cutoff for enrichment
# jmod.long <- "Liver_SV129,Liver_BmalKO"
jmod.long <- "Kidney_SV129"
jmod.long <- "Kidney_SV129,Kidney_BmalKO"

comp <- 1
genes <- as.character(subset(fits.long.filt, model %in% jmod.long)$gene)
s <- SvdOnComplex(subset(dat.freq, gene %in% genes), value.var = "exprs.transformed")
eigens <- GetEigens(s, period = 24, comp = comp, label.n = 25, eigenval = TRUE, adj.mag = TRUE, constant.amp = 4, peak.to.trough = TRUE)
print(eigens$u.plot)


# Highlight genes assigned to GO terms ------------------------------------

genes.bg <- as.character(subset(fits.long.filt)$gene)
genes.fg <- as.character(subset(fits.long.filt, model %in% jmod.long)$gene)

enrichment <- GetGOEnrichment(genes.bg, genes.fg, fdr.cutoff = jcutff, ontology = jonto, show.top.n = 10)

m1 <- ggplot(subset(enrichment, !is.na(Term)), aes(x = Term, y = minuslogpval)) + geom_bar(stat = "identity") + 
  ylab("-log10(P-value), Fisher's exact test") + 
  xlab("") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle(jmod.long)
print(m1)

# annotate dat.freq with GO terms for each gene (take most significant enrichment)

# label clocks to based on GOterm
go.hash <- hash()
for (i in seq(nrow(enrichment))){
  genes.vec <- enrichment$genes[[i]]
  term <- as.character(enrichment$Term[[i]])
  for (g in genes.vec){
    if (is.null(go.hash[[g]])){
      go.hash[[g]] <- term
    } 
  }
}

dat.sub <- subset(dat.freq, gene %in% genes.fg)
s.sub <- SvdOnComplex(dat.sub, value.var = "exprs.transformed")

# label dat.freq with goterm
dat.sub$term <- sapply(as.character(dat.sub$gene), function(g){
  if (!is.null(go.hash[[g]])){
    return(go.hash[[g]])
  } else {
    return(NA)
  }
})

eig <- GetEigens(s.sub, period = 24, comp = comp, label.n = 25, eigenval = TRUE, adj.mag = TRUE, constant.amp = 4, peak.to.trough = TRUE)

omega <- 2 * pi / 24
ampscale <- 2
vec.complex <- eig$eigensamp 
labels <- names(vec.complex)

dat <- data.frame(amp = Mod(vec.complex) * ampscale,
                 phase = ConvertArgToPhase(Arg(vec.complex), omega = omega),
                 label = labels)

# add term BEFORE filtering out labels
dat$term <- as.factor(sapply(as.character(dat$label), function(g){
  if (g == ""){
    return("")
  }
  if (!is.null(go.hash[[g]])){
    return(go.hash[[g]])
  } else {
    return("")
  }
}))
# make "" the first term for levels (so you get colours of dots wihtout the labels 
if (any(dat$term == "")){
  jterms <- as.character(enrichment$Term)
  dat$term <- factor(as.character(dat$term), levels = c("", jterms))
}

top.hits <- 25
top.amps <- as.character(head(dat[order(dat$amp, decreasing = TRUE), ], n = top.hits)$label)
dat$label <- sapply(as.character(dat$label), function(l) ifelse(l %in% top.amps, yes = l, no = ""))
# label only top genes




amp.max <- ceiling(max(dat$amp) * 2) / 2
if (amp.max <= 1){
  amp.step <- 0.5
} else {
  amp.step <- 1
}

xlab <- ""
ylab <- "Log2 Fold Change"
constant.amp <- 6
jtitle <- ""

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbPalette <- cbPalette[1:length(levels(dat$term))]
# change "" term to light gray
cbPalette[levels(dat$term) == ""] <- "gray85"

m <- ggplot(data = dat, aes(x = amp, y = phase, label = label, colour = term)) + 
  geom_point(size = 1) +
  coord_polar(theta = "y") + 
  xlab(xlab) +
  ylab(ylab) +
  ggtitle(jtitle) +
  scale_y_continuous(limits = c(0, 24), breaks = seq(6, 24, 6)) + 
  scale_x_continuous(limits = c(0, amp.max), breaks = seq(0, amp.max, length.out = 2)) + 
  theme_bw() + 
  geom_vline(xintercept = seq(0, amp.max, length.out = 2), colour = "grey50", size = 0.2, linetype = "dashed") +
  geom_hline(yintercept = seq(6, 24, by = 6), colour = "grey50", size = 0.2, linetype = "solid") +
  theme(panel.grid.major = element_line(size = 0.5, colour = "grey"), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position="bottom",
        panel.border = element_blank(),
        legend.key = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank())
# expand_limits(x = 0)  # we want center to be 0

# add text
df.txt <- subset(dat, label != "")
if (constant.amp != FALSE){
  # m <- m + geom_text_repel(data = df.txt, aes(x = amp, y = phase, label = label), size = constant.amp, colour = "black")
  m <- m + geom_text_repel(data = df.txt, aes(x = amp, y = phase, label = label), size = constant.amp)
} else {
  # m <- m + geom_text_repel(data = df.txt, aes(x = amp, y = phase, size = amp, label = label), colour = "black")
  m <- m + geom_text_repel(data = df.txt, aes(x = amp, y = phase, size = amp, label = label))
}
print(m + scale_colour_manual(values=cbPalette) + ggtitle(jonto))

# PlotEnrichmentGenes(dat.freq, enrichment, max.genes = 40, row.i = "max")

