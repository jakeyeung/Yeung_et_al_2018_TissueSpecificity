# cluster_tissues_by_TFs.R
# Jake Yeung
# March 9 2015

library(ggplot2)
# library(plyr)
library(dplyr)
library(reshape2)
library(gplots)

# Paths -------------------------------------------------------------------

# define dirs
data.dir <- "data"
tf.fname <- "motifs_and_TFs.list"
tf.path <- file.path(data.dir, tf.fname)

# Functions ---------------------------------------------------------------

source("scripts/functions/LoadArrayRnaSeq.R")
source("scripts/functions/GetTFs.R")
source("scripts/functions/FirstColToRowname.R")
source("scripts/functions/ClusterTFsFunctions.R")
source("scripts/functions/PlotGeneAcrossTissues.R")


# Set pallette for heatmap colors -----------------------------------------


# creates a own color palette from red to green
# my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)
my.palette <- colorRampPalette(c("black", "yellow"))(n = 299)

# # (optional) defines the color breaks manually for a "skewed" color transition
col.breaks = c(seq(0, 5, length=150),  # black
               seq(5, 20, length=151))  # yellow


# Main --------------------------------------------------------------------

dat <- LoadArrayRnaSeq()

# plot distribution of expression 
ggplot(aes(x = exprs, fill = experiment), data = dat) +
  geom_density(alpha = 0.3)

tfs <- GetTFs()

dat.tfs <- subset(dat, gene %in% tfs & experiment == "rnaseq")

dat.tfs.rnaseq <- subset(dat.tfs, experiment="rnaseq")

dat.tfs.max <- dat.tfs.rnaseq %>%
  group_by(gene) %>%
  do(MaxExprsAcrossTissuesDplyr(.))

dat.tfs.mat <- dcast(dat.tfs.max, gene ~ tissue, value.var = "max.exprs")

dat.tfs.mat <- as.matrix(FirstColToRowname(dat.tfs.mat))

# binarize
cutoff <- 5
dat.tfs.mat.binary <- Binarize(dat.tfs.mat, cutoff)
dat.tfs.mat.binary <- FilterAllOnesOrZeros(dat.tfs.mat.binary)


# heatmap that shiet
heatmap.2(dat.tfs.mat.binary,
          # cellnote = dat.tfs.mat,  # same data set for cell labels
          main = paste0("exprs cutoff=", cutoff), # heat map title
          # notecol="black",      # change font color of cell labels to black
          density.info='density',  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          # margins =c(5, 45),     # widens margins around plot
          dendrogram = "both",
          cexRow = 0.3,
          col=my.palette,       # use on color palette defined earlier 
          # breaks=col.breaks    # enable color transition at specified limits
          # dendrogram="col",     # only draw a row dendrogram
          # Colv="NA")            # turn off column clustering
)
# heatmap that shiet
heatmap.2(dat.tfs.mat,
          # cellnote = dat.tfs.mat,  # same data set for cell labels
          main = "TF exprs", # heat map title
          # notecol="black",      # change font color of cell labels to black
          density.info='density',  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(5, 45),     # widens margins around plot
          dendrogram = "both",
          cexRow = 0.1,
          col=my.palette,       # use on color palette defined earlier 
          # breaks=col.breaks    # enable color transition at specified limits
          # dendrogram="col",     # only draw a row dendrogram
          # Colv="NA")            # turn off column clustering
)

# for fun: Plot Liver vs. Adrenal gland TF expression
dat.liver_v_adr.binary <- dat.tfs.mat.binary[, c("Liver", "Adr")]
dat.liver_v_adr.binary <- FilterAllOnesOrZeros(dat.liver_v_adr.binary)
heatmap.2(dat.liver_v_adr.binary,
          # cellnote = dat.tfs.mat,  # same data set for cell labels
          main = paste0("exprs cutoff=", cutoff), # heat map title
          # notecol="black",      # change font color of cell labels to black
          density.info='density',  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(7, 20),     # widens margins around plot
          dendrogram = "both",
          cexRow = 1,
          col=my.palette,       # use on color palette defined earlier 
          # breaks=col.breaks    # enable color transition at specified limits
          # dendrogram="col",     # only draw a row dendrogram
          # Colv="NA")            # turn off column clustering
)


# Which ones are specific to only ONE tissues -----------------------------

sums <- rowSums(dat.tfs.mat.binary)

outfile <- paste0("plots/tissue_specific_rhythmic_genes/tissue_specific_factors.pdf")
pdf(outfile)
for (N in seq(12)){
  tissue.specific <- names(sums[which(sums == N)])
  if (length(tissue.specific) > 0){
    for (g in tissue.specific){
      jtitle = paste(g, N)
      print(PlotGeneAcrossTissues(subset(dat.tfs.rnaseq, gene == g), jtitle = jtitle))
    }
  }
}
dev.off()
