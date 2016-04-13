# Jake Yeung
# 2016-04-13
# assign_mara_promoters_with_gc_at

library(hash)


# Load --------------------------------------------------------------------


N.path <- "/home/yeung/projects/tissue-specificity/data/sitecounts/motevo/sitecount_matrix"
N.path.meta <- "/home/yeung/projects/tissue-specificity/data/sitecounts/motevo/mara_promoters_gene_name_association.bed"

N <- read.table(N.path, header=T)
N.meta <- read.table(N.path.meta, 
                     col.names = c("chromo", "start", "end", 
                                   "motevo.id", "length", "strand", "gene", 
                                   "transcript.type", "gc.at", "promoter.id"))

common.promoters <- intersect(as.character(N.meta$promoter.id), as.character(N$Promoter))

N <- subset(N, Promoter %in% common.promoters)

c# Hash and rewrite --------------------------------------------------------

gc.hash <- hash(as.character(N.meta$promoter.id), as.character(N.meta$gc.at))

gc.at.vec <- sapply(as.character(N$Promoter), function(p){
  gc.at <- gc.hash[[p]]
  if (is.null(gc.at)) return(NA)
  return(gc.at)
})

length(which(is.na(gc.at.vec)))  # check


# Rewrite -----------------------------------------------------------------


n.motifs <- ncol(subset(N, select=-Promoter))
gc.vec <- rep(c(0, 1), each=n.motifs)
at.vec <- rep(c(1, 0), each=n.motifs)
gcat.lst <- list("at"=at.vec, "gc"=gc.vec)

N.atgc <- cbind(subset(N, select=-Promoter), subset(N, select=-Promoter))
cname.suff <- rep(c("at", "gc"), each=n.motifs)

colnames(N.atgc) <- paste(colnames(N.atgc), cname.suff, sep="_")

N.atgc <- as.matrix(N.atgc)

for (i in seq(nrow(N.atgc))){
  vec <- gcat.lst[[gc.at.vec[i]]]
  N.atgc[i, ] <- N.atgc[i, ] * vec
}

# add Promoter ID back
# N.atgc <- as.data.frame(cbind(N$Promoter, N.atgc))

# Write to table ----------------------------------------------------------

# convert PromoterID to gene
gene.hash <- hash(as.character(N.meta$promoter.id), as.character(N.meta$gene))

gene <- sapply(as.character(N$Promoter), function(p) gene.hash[[p]])

N.atgc <- as.data.frame(cbind(gene, N.atgc))

write.table(N.atgc, file = "/home/yeung/projects/tissue-specificity/data/sitecounts/motevo/sitecount_matrix_geneids.atgc", 
            row.names = FALSE, sep="\t", quote=FALSE)


# Run mara ----------------------------------------------------------------

# for tissuewide
# for all genes

# Analyze mara ------------------------------------------------------------

source("scripts/functions/LoadActivitiesPlotSvd.R")
source("scripts/functions/LoadActivitiesLong.R")

act.long <- LoadActivitiesLong("/home/yeung/projects/tissue-specificity/results/MARA/all_genes_atgc/expressed_genes_deseq_int.centeredTRUE")
LoadActivitiesPlotSvd("/home/yeung/projects/tissue-specificity/results/MARA/all_genes_atgc/expressed_genes_deseq_int.centeredTRUE", 1)
PlotActivitiesWithSE(subset(act.long, gene=="RORA.p2_gc"))
PlotActivitiesWithSE(subset(act.long, gene=="RORA.p2_at"))
PlotActivitiesWithSE(subset(act.long, gene=="bHLH_family.p2_gc"))
PlotActivitiesWithSE(subset(act.long, gene=="bHLH_family.p2_at"))

act.long <- LoadActivitiesLong("/home/yeung/projects/tissue-specificity/results/MARA/bic_modules_at_gc/TissueWide.centeredTRUE")
LoadActivitiesPlotSvd("/home/yeung/projects/tissue-specificity/results/MARA/bic_modules_at_gc/TissueWide.centeredTRUE", 1)

PlotActivitiesWithSE(subset(act.long, gene=="HIC1.p2_at"))


LoadActivitiesPlotSvd("/home/yeung/projects/tissue-specificity/results/MARA/bic_modules_at_gc/TissueWide.centeredTRUE", 1)
act.long <- LoadActivitiesLong("/home/yeung/projects/tissue-specificity/results/MARA/bic_modules_at_gc/TissueWide.centeredTRUE")
PlotActivitiesWithSE(subset(act.long, gene=="RORA.p2_gc"))
PlotActivitiesWithSE(subset(act.long, gene=="RORA.p2_at"))
LoadActivitiesPlotSvd("/home/yeung/projects/tissue-specificity/results/MARA/all_genes_atgc/expressed_genes_deseq_int.centeredTRUE", 3)


# Are AT and GC differentially expressed? ---------------------------------

dat.mean <- subset(dat.long, gene %in% as.character(unique(N.meta$gene)) & experiment=="rnaseq") %>%
  group_by(gene, tissue) %>%
  summarise(exprs=mean(exprs))

N.meta.sum <- subset(N.meta) %>%
  group_by(gene) %>%
  summarise(gc.at=paste(unique(sort(gc.at)), collapse = "_"))

gc.hash.gene <- hash(as.character(N.meta.sum$gene), as.character(N.meta.sum$gc.at))

dat.mean$gc.at <- sapply(as.character(dat.mean$gene), function(g) gc.hash.gene[[g]])

ggplot(dat.mean, aes(x = gc.at, y = exprs)) + geom_boxplot() + facet_wrap(~tissue)

ggplot(dat.mean, aes(x = exprs)) + geom_histogram(bins = 100) + facet_wrap(~gc.at)


# Do MARA on AT and GC separately (different lambda) ----------------------

N.gc <- N[which(gc.at.vec=="gc"), ]
N.gc.genes <- gene[which(gc.at.vec=="gc")]
N.at <- N[which(gc.at.vec=="at"), ]
N.at.genes <- gene[which(gc.at.vec=="at")]

write.table(cbind(Gene.ID=N.gc.genes, subset(N.gc, select=-Promoter)), file = "/home/yeung/projects/tissue-specificity/data/sitecounts/motevo/sitecount_matrix_geneids.gc_only", quote=FALSE, sep="\t", row.names=FALSE)
write.table(cbind(Gene.ID=N.at.genes, subset(N.at, select=-Promoter)), file = "/home/yeung/projects/tissue-specificity/data/sitecounts/motevo/sitecount_matrix_geneids.at_only", quote=FALSE, sep="\t", row.names=FALSE)


# Analyze MARA output -----------------------------------------------------

LoadActivitiesPlotSvd("/home/yeung/projects/tissue-specificity/results/MARA/all_genes_at_only/expressed_genes_deseq_int.centeredTRUE", 1)
LoadActivitiesPlotSvd("/home/yeung/projects/tissue-specificity/results/MARA/all_genes_gc_only/expressed_genes_deseq_int.centeredTRUE", 1)

act.long.gc <- LoadActivitiesLong("/home/yeung/projects/tissue-specificity/results/MARA/all_genes_gc_only/expressed_genes_deseq_int.centeredTRUE")
act.long.at <- LoadActivitiesLong("/home/yeung/projects/tissue-specificity/results/MARA/all_genes_at_only/expressed_genes_deseq_int.centeredTRUE")

PlotActivitiesWithSE(subset(act.long.gc, gene=="HNF4A_NR2F1.2.p2"))
PlotActivitiesWithSE(subset(act.long.at, gene=="HNF4A_NR2F1.2.p2"))
PlotActivitiesWithSE(subset(act.long.gc, gene=="RORA.p2"))
PlotActivitiesWithSE(subset(act.long.at, gene=="RORA.p2"))
PlotActivitiesWithSE(subset(act.long.gc, gene=="HIC1.p2"))
PlotActivitiesWithSE(subset(act.long.at, gene=="HIC1.p2"))
