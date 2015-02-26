# FilterBedFiles.R
#
# Given list of genes, filter bed file so it only contains gene of interest.


# Functions ---------------------------------------------------------------

source(file.path("scripts", "functions", "FilterBedFile.R"))
source(file.path("scripts", "functions", "ReadListToVector.R"))
source(file.path("scripts", "functions", "RemoveExtension.R"))


fnames.list <- "~/projects/tissue-specificity/plots/nconds/7_conds_filtered_02_bicw/files.txt"
fnames <- ReadListToVector(fnames.list)
outdir <- "~/projects/tissue-specificity/bedfiles/filtered_02_bicw"

for (fname in fnames){
  bedfile.out <- file.path(outdir, basename(fname))
  bedfile.out <- RemoveExtension(bedfile.out)
  bedfile.out <- paste0(bedfile.out, '.bed')
  FilterBedFile(fname, bedfile.out)
}

# do similar with filtered_genes.txt
bedfile.out <- file.path(outdir, "filtered_genes.bed")
FilterBedFile(file.path(outdir, "filtered_genes.txt", bedfile.out))




