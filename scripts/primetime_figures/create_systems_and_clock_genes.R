# Jake Yeung
# Date of Creation: 2019-02-21
# File: ~/projects/tissue-specificity/scripts/primetime_figures/create_systems_and_clock_genes.R
# Create systems and clock for Colas and Felix

rm(list=ls())

# Load primetime objects --------------------------------------------------

inf <- "/data/shared/jake_data/tissue_specificity/PostPhDFiles/GR_2018_Primetime_Objects.Rdata"
load(inf, v=T)

# Find tissue wide genes --------------------------------------------------

genes.tw <- as.character(subset(fits.long, n.rhyth >= 8)$gene)
genes.tw.wtko <- as.character(subset(fits.long.filt, model %in% c("Liver_SV129,Kidney_SV129"))$gene)


# Assign as clock or system -----------------------------------------------

fits.sub <- subset(fits.long.filt, gene %in% genes.tw & model != "" & amp.avg > 0) 

fits.sum <- fits.sub %>% 
  group_by(model) %>% 
  summarise(count = length(gene)) %>%
  arrange(desc(count))


# Assign to clock or system -----------------------------------------------


clock.or.system <- sapply(fits.sum$model, function(m){
  # if model is Liver_SV129,Kidney_SV129;Liver_BmalKO,Kidney_BmalKO, consider
  # it clock
  if (m == "Liver_SV129,Kidney_SV129;Liver_BmalKO,Kidney_BmalKO"){
    return("clock")
  }
  
  if (grepl("KO", m)){
    return("system")
  } else if (grepl("SV129", m)){
    return("clock")
  } else {
    warning("Neither KO or SV129")
    return(NA)
  }
})
clock.sys.hash <- hash(as.character(fits.sum$model), clock.or.system)

fits.sub$clksys <- sapply(as.character(fits.sub$model), function(m) clock.sys.hash[[m]])

clock.sys.gene.hash <- hash(as.character(fits.sub$gene), fits.sub$clksys)


# Replot tissue-wide, colored by clock or systems

comp <- 1
col.hash.gene <- hash()
for (g in genes.tw){
  col.hash.gene[[g]] <- ifelse(clock.sys.gene.hash[[g]] == "clock", "red", "blue")
}
s.tw <- SvdOnComplex(subset(dat.complex, gene %in% genes.tw), value.var = "exprs.transformed")

eigens.tw <- GetEigens(s.tw, period = 24, comp = comp, adj.mag = TRUE, 
                       eigenval = TRUE,
                       constant.amp = 5, 
                       label.n = Inf, jtitle = "", 
                       peak.to.trough = TRUE, 
                       dot.col = col.hash.gene, 
                       dotsize = 1.5, 
                       dotshape = 18,
                       disable.text = TRUE, 
                       add.arrow = TRUE,
                       disable.repel = TRUE,
                       half.life = 0)

jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)
print(eigens.tw$u.plot)
print(eigens.tw$v.plot)
multiplot(eigens.tw$u.plot, eigens.tw$v.plot, layout = jlayout)

eigens.tw <- GetEigens(s.tw, period = 24, comp = comp, adj.mag = TRUE, 
                       eigenval = TRUE,
                       constant.amp = 5, 
                       label.n = Inf, jtitle = "", 
                       peak.to.trough = TRUE, 
                       dot.col = col.hash.gene, 
                       dotsize = 3, 
                       dotshape = 18,
                       disable.text = TRUE, 
                       add.arrow = TRUE,
                       disable.repel = TRUE,
                       half.life = 0)

jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)
print(eigens.tw$u.plot)
print(eigens.tw$v.plot)
multiplot(eigens.tw$u.plot, eigens.tw$v.plot, layout = jlayout)


# Write to output ---------------------------------------------------------


data.table::fwrite(fits.sub, file = "/data/shared/jake_data/tissue_specificity/PostPhDFiles/tissue_wide_genes_with_clcksys_label.txt", sep = "\t")



