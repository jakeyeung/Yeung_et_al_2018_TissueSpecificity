# explore_site_counts.R
# Saeed has a new sitecounts, let's see what's new.


# Functions ---------------------------------------------------------------

PlotSitecounts <- function(vec, filter.top = 40, title = 'MyTitle'){
  # Given vector of sitecounts for a gene and labels for its motifs,
  # plot a labelled bar graph
  #
  # Remove zeros
  vec <- vec[which(vec > 0)]
  i.ordered <- order(vec, decreasing = TRUE)
  # Filter top 20
  i.ordered <- head(i.ordered, n = filter.top)
  vec.ord <- vec[i.ordered]
  heights <- rev(vec.ord)  # first bar at top
  barplot(height = heights, names.arg = names(heights), main = title, horiz = TRUE, las = 1)
}

# Main --------------------------------------------------------------------


# Load data ---------------------------------------------------------------

N.path <- "data/sitecounts/motevo/sitecount_matrix_geneids"
N.promoterpath <- "data/sitecounts/motevo/mara_promoters_gene_name_association.bed"

N <- read.table(N.path, header = TRUE, sep = '\t')
rnames <- make.names(N$Gene.ID, unique = TRUE)
rownames(N) <- rnames
N$Gene.ID <- NULL

N.promoter <- read.table(N.promoterpath, header = FALSE, sep = '\t')
N.promoter <- data.frame(Gene.ID=rnames, maraid=N.promoter$V4, saeedid=N.promoter$V10)
rownames(N.promoter) <- rnames

# Plot ddc, our favourite gene --------------------------------------------

gene.hit <- "Ddc"  # regulated by ROR
gene.hit <- "Sept9"  # regulated by PPARG MEF2D, TEAD, noncore clocks
gene.hit <- "Steap3"  # Steap3

# par(mar=c(5.1,4.1,4.1,2.1))  # default
par(mar=c(5.1, 12, 4.1, 2.1))

for (gene in rnames[grep(gene.hit, rnames)]){
  promoter <- N.promoter[gene, ]$maraid
  title <- paste(gene, promoter)
  print(title)
  PlotSitecounts(unlist(N[gene, ]), filter.top = 30, title = title)
  # promoter <- N.promoter[[gene]]
  # plot(unlist(N[gene, ]), col = "white", main = paste(c(gene, promoter)))
  # text(unlist(N[gene, ]), labels = names(unlist(N[gene, ])), cex = unlist(N[gene, ]) + 0.01)
}

