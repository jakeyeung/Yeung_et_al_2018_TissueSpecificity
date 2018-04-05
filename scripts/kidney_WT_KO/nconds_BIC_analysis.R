# Jake Yeung
# Date of Creation: 2018-01-20
# File: ~/projects/tissue-specificity/scripts/kidney_WT_KO/nconds_BIC_analysis.R
# Look at BIC genes because the g sweep is bizarre with a single rhythm

library(dplyr)
library(ggplot2)
library(reshape2)

# Load --------------------------------------------------------------------

jmeth <- "g=1001"
jmeth <- "BIC"
inf <- "/home/yeung/projects/tissue-specificity/Robjs/liver_kidney_atger_nestle/nconds_kidney_WTKO/summary/fits.kidney_WTKO.multimethod.long.filtbest.staggeredtimepts.bugfixed.Robj"
load(inf, v=T)

fits <- subset(fits.long.filt, method == jmeth)

load("/home/yeung/projects/tissue-specificity/Robjs/liver_kidney_atger_nestle/nconds_kidney_WTKO/summary/dat.freq.kidneyWTKO.bugfixed.Robj", v=T)

# load RNASeq
load("Robjs/liver_kidney_atger_nestle/dat.long.liverkidneyWTKO.bugfixed.Robj", v=T)
dat.long <- CollapseTissueGeno(dat.long)
removesamps <- TRUE
stagger <- TRUE
if (removesamps){
  if (stagger){
    # Stagger removes ZT2 from Kidney, which is an outlier when we use Kallisto
    dat.long <- StaggeredTimepointsLivKid(dat.long)
  } else {
    dat.long <- SameTimepointsLivKid(dat.long)
  }
}
dat.long <- subset(dat.long, tissue == "Kidney_SV129" | tissue == "Kidney_BmalKO")
dat.long$tissue <- factor(as.character(dat.long$tissue), levels = c("Kidney_SV129", "Kidney_BmalKO"))

dat.mat <- dcast(dat.long, formula = "gene ~ tissue + geno + time", value.var = "exprs")

# add models to dat.mat
# parse out params.list into columns
params.full.lst <- lapply(fits$param.list, function(lst){
  # expect names: "tissueKidney_SV129, tissueKidney_BmalKO, Kidney_SV129.amp, Kidney_BmalKO.amp, Kidney_SV129.phase, Kidney_BmalKO.phase"
  cnames <- c("tissueKidney_SV129", "tissueKidney_BmalKO", "Kidney_SV129.amp", "Kidney_BmalKO.amp", "Kidney_SV129.phase", "Kidney_BmalKO.phase", 
              "Kidney_SV129,Kidney_BmalKO.amp", "Kidney_SV129,Kidney_BmalKO.phase")
  param.vals <- sapply(cnames, function(x){
    if (x %in% names(lst)){
      return(lst[[x]])
    } else {
      return(NA)
    }
  })
})
params.full <- as.data.frame(do.call(rbind, params.full.lst))
params.full$gene <- as.character(fits$gene)
# rename cnames to be more interpretable
colnames(params.full) <- gsub(pattern = "tissue", replacement = "intercept", colnames(params.full))
dat.mat.params <- dplyr::inner_join(dat.mat, subset(fits, select = c(gene, model, weight)), by = "gene")
dat.mat.params <- dplyr::inner_join(dat.mat.params, params.full, by = "gene")

head(dat.mat.params)

# How many models?  -------------------------------------------------------

fits.sum <- fits %>%
  group_by(model) %>%
  summarise(number.of.genes = length(gene))

m.bar <- ggplot(fits.sum, aes(x = model, y=number.of.genes)) + geom_bar(stat = "identity") + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

m.bar.noflat <- ggplot(subset(fits.sum, model != ""), aes(x = model, y=number.of.genes)) + geom_bar(stat = "identity") + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Plot SVD ----------------------------------------------------------------

i <- 1
jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)

genes.clock <- as.character(subset(fits, model %in% c("Kidney_SV129"))$gene)
s <- SvdOnComplex(subset(dat.freq, gene %in% genes.clock), value.var = "exprs.transformed")
eigens.clock <- GetEigens(s, period = 24, comp = i, label.n = 30, eigenval = TRUE, adj.mag = TRUE, constant.amp = 4, peak.to.trough = TRUE)


genes.system <- as.character(subset(fits, model %in% c("Kidney_SV129,Kidney_BmalKO"))$gene)
s <- SvdOnComplex(subset(dat.freq, gene %in% genes.system), value.var = "exprs.transformed")
eigens.system <- GetEigens(s, period = 24, comp = i, label.n = 30, eigenval = TRUE, adj.mag = TRUE, constant.amp = 4, peak.to.trough = TRUE)

genes.full <- as.character(subset(fits, model %in% c("Kidney_SV129;Kidney_BmalKO"))$gene)
s <- SvdOnComplex(subset(dat.freq, gene %in% genes.full), value.var = "exprs.transformed")
eigens.full <- GetEigens(s, period = 24, comp = 1, label.n = 30, eigenval = TRUE, adj.mag = TRUE, constant.amp = 4, peak.to.trough = TRUE)

genes.KO <- as.character(subset(fits, model %in% c("Kidney_BmalKO"))$gene)
s <- SvdOnComplex(subset(dat.freq, gene %in% genes.KO), value.var = "exprs.transformed")
eigens.KO <- GetEigens(s, period = 24, comp = 1, label.n = 30, eigenval = TRUE, adj.mag = TRUE, constant.amp = 4, peak.to.trough = TRUE)

# plot
outdir <- "/home/yeung/projects/tissue-specificity/plots/kidney_WTKO"
pdf(file.path(outdir, paste0("rhythmic_modules_svd_method.", jmeth, ".pdf")), useDingbats = FALSE)
  print(m.bar)
  print(m.bar.noflat)
  eigens.clock$u.plot + ggtitle("Clock-regulated genes") + theme(plot.title = element_text(size = 12))
  eigens.clock$v.plot + ggtitle("Clock-regulated genes") + theme(plot.title = element_text(size = 12))
  
  eigens.system$u.plot + ggtitle("System-driven genes") + theme(plot.title = element_text(size = 12))
  eigens.system$v.plot + ggtitle("System-driven genes") + theme(plot.title = element_text(size = 12))
  
  eigens.full$u.plot + ggtitle("Full model") + theme(plot.title = element_text(size = 12))
  eigens.full$v.plot + ggtitle("Full model") + theme(plot.title = element_text(size = 12))
  
  eigens.KO$u.plot + ggtitle("Bmal1KO-specific rhythms") + theme(plot.title = element_text(size = 12))
  eigens.KO$v.plot + ggtitle("Bmal1KO-specific rhythms") + theme(plot.title = element_text(size = 12))
dev.off()


# write table

write.table(dat.mat.params, file = file.path(outdir, paste0("rnaseq_kidney_log_tpm_and_models.method.", jmeth, ".txt")), quote = FALSE, sep = "\t", row.names = FALSE)

