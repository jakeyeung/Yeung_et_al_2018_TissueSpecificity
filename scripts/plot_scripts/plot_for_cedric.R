# 2016-03-29
# Cedric requests a plot of genes

# copied from email 
jgenes_str <- "SOCS3
FOS
PDK4
ANKRD1
ANGPTL4
EGR1
THBS1
ZFP36
MYC
RGS16
CXCL10
CDKN1A
GADD45A
PFKFB3
ATF3
LGALS17A
C11orf96
RP11-613D13.4
DUSP2"

jgenes.upper <- strsplit(jgenes_str, "\n")[[1]]

# Make lower case on non-first letter
jgenes <- sapply(jgenes.upper, function(jgene){
  paste(toupper(substr(jgene, 1, 1)), tolower(substr(jgene, 2, nchar(jgene))), sep="")
})

load("Robjs/dat.long.fixed_rik_genes.Robj", v=T)

dat.sub <- subset(dat.long, gene %in% jgenes)
source("scripts/functions/PlotGeneAcrossTissues.R"
pdf("plots/gene_expressions/cedric_genes.pdf")
for (jgene in jgenes){
  dat.gene <- subset(dat.sub, gene == jgene)
  if (nrow(dat.gene) > 0){
    print(PlotGeneAcrossTissues(subset(dat.sub, gene == jgene)))
  } else {
    print(paste("Skipping", jgene))
  }
}
dev.off()
