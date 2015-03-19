# load_sobel_wang_data.R
# Jake Yeung
# Date: 2015-03-18
# from http://jefworks.com/gene-id-conversion-in-r-using-biomart/


library("biomaRt")

Transcript2Gene <- function(gene.list) {
  mart.obj <- useMart(biomart = 'ensembl', dataset = 'mmusculus_gene_ensembl')
  gos <- getBM(gene.list,attributes=c("ensembl_transcript_id", "external_gene_name"),
               filters=c("ensembl_transcript_id"),
               mart=mart.obj)
  gl <- gos[match(gene.list, gos[,1]), 2]
  ## if not found, then keep ENSEMBL name
  print(paste0("Could not match ", length(gl[is.na(gl)]), " genes."))
  # gl[is.na(gl)] <- gene.list[is.na(gl)]
  return(gl)
}

AppendGeneID <- function(dat){
  genes.tens <- rownames(dat)
  Gene.ID <- Transcript2Gene(genes.tens)
  dat <- cbind(Gene.ID, dat)
  return(dat)
}

# Load RData --------------------------------------------------------------

sobel.obj.path <- "data/from_labmates/matrix_for_jake.RData"  # pm_WT and x33_WT
wang.obj.path <- "data/from_labmates//Elastic_TFs_Input_Pol2_CM_library_final.Rdata"  # mot_pol2 and res.pol2.sel
dhs.meta.path <- "site_count_matrices//annotated_DNAse_peaks.txt"
load(sobel.obj.path)
load(wang.obj.path)

dhs.annots <- read.table(dhs.meta.path, header = TRUE)

coords <- rownames(pm_WT)
genes <- subset(dhs.annots, index %in% coords, select = featurename)
coords_to_genes <- data.frame(coord = coords, featurename = genes)

genes.uniq <- unique(unlist(genes))  # transcript ID actually
pm_WT.collapsed <- matrix(data = NA, nrow = length(genes.uniq), ncol = ncol(pm_WT), dimnames = list(genes.uniq, colnames(pm_WT)))

for (gene in genes.uniq){
  coords <- coords_to_genes[which(coords_to_genes$featurename == gene), 1]
  counts <- pm_WT[coords, ]
  counts.collapsed <- apply(counts, 2, sum)
  pm_WT.collapsed[gene, ] <- counts.collapsed
}

pm_WT.genes <- Transcript2Genes(genes.uniq)  # contains some NAs
pm_WT.collapsed.genenames <- as.data.frame(AppendGeneID(pm_WT.collapsed))

write.table(pm_WT.collapsed.genenames,
            file = "/home/yeung/projects/tissue-specificity/site_count_matrices/GLM/Wang/sitecount_matrix.full",
            quote = FALSE,
            sep = "\t",
            row.names = FALSE)
save(pm_WT.collapsed, file = "/home/yeung/projects/tissue-specificity/site_count_matrices/GLM/Wang/pm_WT.Robj")

# # Collapse same names to a single name ------------------------------------
# 
# mot_pol2 <- as.data.frame(AppendGeneID(mot_pol2))
# res.pol2.sel <- as.data.frame(AppendGeneID(res.pol2.sel))
# 
# # make rownames
# rnames <- make.names(res.pol2.sel$Gene.ID, unique = TRUE)
# rownames(res.pol2.sel) <- rname
# res.pol2.sel$phase <- as.numeric(levels(res.pol2.sel$phase))[as.integer(res.pol2.sel$phase)]
# 
# # save Jingkui's output to file
# write.table(mot_pol2, 
#             file = "/home/yeung/projects/tissue-specificity/site_count_matrices/GLM/Wang/sitecount_matrix",
#             quote = FALSE,
#             sep = "\t",
#             row.names = FALSE)
