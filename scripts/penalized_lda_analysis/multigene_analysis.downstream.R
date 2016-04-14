# 2016-04-13
# Jake Yeung
# multigene_analysis.downstream.R

setwd("/home/yeung/projects/tissue-specificity")
  
library(wordcloud)

# Load --------------------------------------------------------------------

source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/PlotUCSC.R")

load("Robjs/penalized_lda_mats.5000.Liver.cutoff2.Robj", v=T)
load("Robjs/dat.long.fixed_rik_genes.Robj", v=T)
load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.15.phase_sd_maxdiff_avg.Robj", v=T)

# Sort by ROR occurances --------------------------------------------------

top.hits <- mat.fg[order(mat.fg$RORA.p2, decreasing = TRUE), c("peak", "gene", "RORA.p2")]
top.hits <- top.hits[which(top.hits$RORA.p2 > 0), ]

jgene <- "Agtpbp1"
jpeak <- "chr13:59481037-59481537"

# jgene <- "Gm10447"
# jpeak <- "chr11:53452634-53453134"
jgene <- "Leap2"
jpeak <- "chr11:53452634-53453134"

jgene <- "Agtpbp1"
jpeak <- "chr13:59481037-59481537"

jgene <- "Fbxo44"
jpeak <- "chr4:148147391-148147891"

outf <- "/home/yeung/projects/tissue-specificity/plots/penalized_lda/multigeen_rora_analysis.5000.cutoff2.ucsc.pdf"
if (!file.exists(outf)){
  pdf(outf)
  for (i in seq(nrow(top.hits))){
    jgene <- as.character(top.hits[i, ]$gene[[1]])
    jpeak <- as.character(top.hits[i, ]$peak[[1]])
    print(paste(jgene, jpeak))
    
    print(PlotGeneAcrossTissues(subset(dat.long, gene == jgene)))
    
    # show counts at peak ranked by motifs
    jsub <- subset(mat.fg, peak == jpeak & gene == jgene, select=c(-peak, -gene))
    jsub <- jsub[which(jsub > 0)]
    jsub <- t(as.matrix(jsub[order(jsub, decreasing=TRUE)]))
    textplot(x = seq(length(jsub)), y = jsub, words = rownames(jsub), main=paste(jgene, jpeak))
    
    # do for all liver-specific peaks near gene
    
    jsub.gene <- colSums(subset(mat.fg, gene == jgene, selec=c(-peak, -gene)))
    jsub.gene <- jsub.gene[which(jsub.gene > 0)]
    jsub.gene <- as.matrix(jsub.gene[order(jsub.gene, decreasing=TRUE)])
    textplot(x = seq(length(jsub.gene)), y = jsub.gene, words=rownames(jsub.gene), main=paste(jgene, "all DHS"))
  }
  dev.off()
} else {
  print(paste("Skipping outf", outf))
}


# Plot UCSC ---------------------------------------------------------------
print("Printing UCSC")
bed <- lapply(as.character(top.hits$peak), CoordToBed); bed <- do.call(rbind, bed)

bedToUCSC(toPlot = bed, outpdf = "plots/penalized_lda/multigeen_rora_analysis.5000.cutoff2.ucsc.pdf", leftwindow = 10000, rightwindow = 10000, 
          theURL = "http://genome-euro.ucsc.edu/cgi-bin/hgTracks?hgsid=213407622_Xnr9Qa6ZR3hPPSc3pFsVhHrz0Etf&hgt.psOutput=on", euro.mirror=TRUE)

