# mara_motif_to_genenames.R
# Jake Yeung
# 2016-03-07

source("scripts/functions/GetTFs.R")
source('scripts/functions/RemoveP2Name.R')
source("scripts/functions/PlotGeneAcrossTissues.R")

load("Robjs/dat.long.fixed_rik_genes.Robj", v=T)


# Load --------------------------------------------------------------------

# from analyze_mara.complex_svd.R filtering for tissuewide module
hits <- c('HIC1','bHLH_family','NFIL3','PAX5','NFIX','TFDP1','ATF2','HSF1.2','NFE2L2',
          'RORA','MYBL2','EP300','HLF','NR3C1','SPZ1','SRY','ATF4','ZBTB6','POU3F1..4','FOXQ1','SRF')

tfs.genes <- GetTFs(split.commas = FALSE)
tfs.motif <- GetTFs(get.motifs = TRUE)
tfs.motif <- sapply(tfs.motif, RemoveP2Name)

tfs.df <- data.frame(motif = tfs.motif, genes = tfs.genes); rownames(tfs.df) <- NULL

hits.genes <- as.character(subset(tfs.df, motif %in% hits)$gene)

hits.genes.flat <- unlist(strsplit(hits.genes, split=","))

pdf("plots/gene_expressions/tissuewide_regulators_transcription.pdf")
for (jgene in hits.genes.flat){
  jsub <- subset(dat.long, gene == jgene)
  if (nrow(jsub) > 0){
    m <- PlotGeneAcrossTissues()
    print(m)
  }
}
dev.off()