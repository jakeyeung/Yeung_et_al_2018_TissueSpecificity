# load_sobel_wang_data.R
# Jake Yeung
# Date: 2015-03-18
# from http://jefworks.com/gene-id-conversion-in-r-using-biomart/

setwd("~/projects/tissue-specificity")

source("scripts/functions//BiomartFunctions.R")

# Load RData --------------------------------------------------------------

sobel.obj.path <- "data/from_labmates/matrix_for_jake.RData"  # pm_WT and x33_WT
wang.obj.path <- "data/from_labmates//Elastic_TFs_Input_Pol2_CM_library_final.Rdata"  # mot_pol2 and res.pol2.sel
dhs.meta.path <- "site_count_matrices//annotated_DNAse_peaks.txt"
load(sobel.obj.path)
load(wang.obj.path)

dhs.annots <- read.table(dhs.meta.path, header = TRUE)

coords <- rownames(pm_WT)
genes.transcript <- subset(dhs.annots, index %in% coords, select = featurename)
genes <- Transcript2Gene(unlist(genes.transcript), return.original=TRUE)
coords_to_genes <- data.frame(coord = coords, featurename = genes)

genes.uniq <- unique(unlist(genes))  # transcript ID converted to gene name
pm_WT.collapsed <- matrix(data = NA, nrow = length(genes.uniq), ncol = ncol(pm_WT), dimnames = list(genes.uniq, colnames(pm_WT)))

for (gene in genes.uniq){
  coords <- coords_to_genes[which(coords_to_genes$featurename == gene), 1]
  counts <- pm_WT[coords, ]
  counts.collapsed <- apply(counts, 2, sum)
  pm_WT.collapsed[gene, ] <- counts.collapsed
}

# pm_WT.genes <- Transcript2Gene(genes.uniq)  # contains some NAs
# pm_WT.collapsed.genenames <- as.data.frame(AppendGeneID(pm_WT.collapsed))

Gene.ID <- rownames(pm_WT.collapsed)
pm_WT.collapsed <- cbind(Gene.ID, pm_WT.collapsed)

write.table(pm_WT.collapsed,
            file = "/home/yeung/projects/tissue-specificity/site_count_matrices/GLM/Wang/sitecount_matrix.full.genenames",
            quote = FALSE,
            sep = "\t",
            row.names = FALSE)
save(pm_WT.collapsed, file = "/home/yeung/projects/tissue-specificity/site_count_matrices/GLM/Wang/pm_WT.Robj")
