# 2015-05-19
# filter_kallisto_for_start_sites.R

library(dplyr)
library(gplots)
library(reshape2)

# Functions ---------------------------------------------------------------

source("scripts/functions/LoadKallisto.R")
source("scripts/functions/AlternativeFirstExonsFunctions.R")
source("scripts/functions/LoadArrayRnaSeq.R")
source("scripts/functions/FitRhythmic.R")
source("scripts/functions/PlotSitecounts.R")
source("scripts/functions/LoadSitecounts.R")

GetOverlaps <- function(starts, ends){
  for (i in 1:length(starts)){
    merged_interval <- c(starts[i], ends[i])
    mygroup <- 
    for (j in (i + 1):length(starts)){
      if (ends[i] > starts[j]){
        # overlap
        merged_interval <- c(starts[i], ends[j])
        
      }
    }
  }
}

MergeOverlappingRegions <- function(dat){
  # Input: dat with columns chromo, start, end.
  # Output: dat with column combined_exon
  
  chromo <- as.character(dat$chromo[1])
  
  if (nrow(dat) == 1){
    # don't do anything if gene only contains one transcript
    dat$merged_intervals <- paste(chromo, paste(dat$start, dat$end, sep = "-"), sep = ":")
    dat$merged_transcripts <- dat$transcript_id
    return(dat)
  }
  
  transcripts <- as.character(dat$transcript_id)
  
  merged_intervals <- rep(NA, nrow(dat))
  merged_transcripts <- rep(NA, nrow(dat))
  
  indices_all = c()  # don't know how big this will big
  for (i in 1:(nrow(dat) - 1)){  # in j we will look at next row, so only iterate until 2nd last row
    if (i %in% indices_all) next
    merged_interval <- c(dat$start[i], dat$end[i])
    merged_transcript <- transcripts[i]
    indices <- i
      for (j in (i + 1):nrow(dat)){
        if (dat$end[i] > dat$start[j]){
          # overlap, lengthen the end of the merged_interval
          merged_interval[2] <- dat$end[j]
          merged_transcript <- c(merged_transcript, transcripts[j])
          indices <- c(indices, j)
        }
      }
    merged_interval <- paste(merged_interval, collapse = "-")
    merged_intervals[indices] <- paste(chromo, merged_interval, sep = ':')
    merged_transcripts[indices] <- paste(merged_transcript, collapse = ",")
    indices_all <- c(indices_all, indices)
  }
  dat$merged_intervals <- merged_intervals
  dat$merged_transcripts <- merged_transcripts
  return(dat)
}

resetPar <- function() {
  dev.new()
  op <- par(no.readonly = TRUE)
  dev.off()
  op
}

# Load data ---------------------------------------------------------------

tpm.merged <- LoadKallisto(path.kallisto = "data/alternative_exon_usage/abundance.merged.annotated.sorted.pythonmerged.bed")
dat.long <- LoadArrayRnaSeq()

# find tissue-specific circadian genes ------------------------------------------------------------

dat.rhyth <- FitRhythmicDatLong(dat.long)
dat.rhyth$is.rhythmic <- mapply(IsRhythmic2, pval = dat.rhyth$pval, avg.exprs = dat.rhyth$int.rnaseq)
dat.rhyth <- dat.rhyth %>%
  group_by(gene) %>%
  do(IsTissueSpecific2(., pval.min = 1e-5, pval.max = 0.05, cutoff = 6))

tissue.spec.circ.genes <- unique(dat.rhyth[which(dat.rhyth$is.tissue.spec.circ == TRUE), ]$gene)

# Calculate isoform usage -------------------------------------------------

tpm.afe <- tpm.merged %>%
  group_by(gene_name, tissue, time) %>%
  mutate(tpm_normalized = CalculateFractionIsoformUsage(tpm, pseudocount = 1))

# Find associations between fractional isoform usage and rhythmici --------

# FIRST: assign genes as "rhythmic" or "not rhythmic"
tpm.afe.filt <- subset(tpm.afe, gene_name %in% tissue.spec.circ.genes)
rhythmic.dic <- setNames(object = dat.rhyth$is.rhythmic, nm = paste(dat.rhyth$gene, dat.rhyth$tissue, sep = ';'))
tpm.afe.filt$is.rhythmic <- rhythmic.dic[paste(tpm.afe.filt$gene_name, tpm.afe.filt$tissue, sep = ';')]
# DONE FIRST

# 2: Associate tissue-specific rhythms with fractional isoform usage
fit.afe <- tpm.afe.filt %>%
  filter(!is.na(is.rhythmic)) %>%  # filter out NA's. Onl consider TRUE and FALSE
  group_by(gene_name, transcript_id) %>%
  do(ModelRhythmicity(., jformula = tpm_normalized ~ is.rhythmic))

# 3. summarize by choosing the top for each gene
fit.afe.summary <- fit.afe %>%
  group_by(gene_name) %>%
  do(SubsetMinPval(jdf = .))

fit.afe.summary$pval.adj <- p.adjust(fit.afe.summary$pval)

# How many make it past threshold?
pval.adj.thres <- 0.05

genes.tested <- unique(fit.afe.summary$gene_name)
n.hits <- length(which(fit.afe.summary$pval.adj <= pval.adj.thres))

sprintf("%s hits found out of %s tissue-specific circadian genes tested. %f percent", 
        n.hits, length(genes.tested), 100 * (n.hits / length(genes.tested)))


# Sanity checks -----------------------------------------------------------

# show top hits (hide transcript_id because it's a huge column)
top.hits <- head(data.frame(fit.afe.summary[order(fit.afe.summary$pval), ]), n = n.hits)
(subset(top.hits, select = -transcript_id))


PlotDiagnostics(dat.tpm, dat.long, "Srrm2", "ENSMUST00000191385,ENSMUST00000190293")
PlotDiagnostics(dat.tpm, dat.long, "Adipoq", "ENSMUST00000023593")
PlotDiagnostics(dat.tpm, dat.long, "Myo1b", "ENSMUST00000144694")
PlotDiagnostics(dat.tpm, dat.long, "Tcp1", "ENSMUST00000089024,ENSMUST00000151287,ENSMUST00000143961,ENSMUST00000129632,ENSMUST00000133003,ENSMUST00000151715")


# start <- Sys.time()
# # 18 minutes: super slow. Comment out if not needed.
# pdf("plots/alternative_exon_usage/kallisto_diagnostics.start_site_merged.top100.pdf", width = 21, height = 7, paper = "USr")
# jgenes <- as.character(top.hits$gene_name); jtranscripts <- as.character(top.hits$transcript_id)
# dat.tpm <- subset(tpm.afe.filt, gene_name %in% jgenes); dat.arrayrnaseq <- subset(dat.long, gene %in% jgenes)
# mapply(PlotDiagnostics, jgene = jgenes, jtranscript = jtranscripts, 
#        MoreArgs = list(dat.tpm = dat.tpm, dat.arrayrnaseq = dat.arrayrnaseq))
# dev.off()
# print(Sys.time() - start)


# What is the distribution of tissues in tissue-specific circadian --------

dat.sub <- subset(dat.rhyth, is.tissue.spec.circ == TRUE & is.rhythmic == TRUE)

tissue_distributions <- dat.sub %>%
  group_by(tissue) %>%
  summarise(count = length(tissue)) %>%
  arrange(desc(count))

# rearrange factor levels
tissue_distributions$tissue <- factor(tissue_distributions$tissue, 
                                      levels = as.character(tissue_distributions$tissue))

ggplot(data = tissue_distributions, aes(x = tissue, y = count)) + 
  geom_bar(stat = "identity") +
  ggtitle("Number of genes containing tissue-specific circadian rhythms")


# What is distribution of hits? -------------------------------------------

dat.sub <- subset(dat.rhyth, gene %in% top.hits$gene_name & is.rhythmic == TRUE)

tissue_distributions <- dat.sub %>%
  group_by(tissue) %>%
  summarise(count = length(tissue)) %>%
  arrange(desc(count))

# rearrange factor levels
tissue_distributions$tissue <- factor(tissue_distributions$tissue, 
                                      levels = as.character(tissue_distributions$tissue))

ggplot(data = tissue_distributions, aes(x = tissue, y = count)) + 
  geom_bar(stat = "identity") +
  ggtitle("Number of genes containing tissue-specific circadian rhythms associated with alternative promoters")


# Explore sitecounts for hits ---------------------------------------------

N.list <- LoadSitecounts(gene_ids=FALSE)  # list of N and N.promoter
N <- as.matrix(N.list$N)
N.annot <- N.list$N.promoter
if (exists(x = "N") & exists(x = "N.promoter")) rm(N.list)

N.long <- LoadSitecountsPromotersLong()

jgene <- "Upp2"
jgene <- "Slc45a3"
jgene <- "Ddc"
jgene <- "Insig2"
jgene <- "Megf10"  # false positive?
jgene <- "Tnfaip8"
jgene <- "Fgf1"
jgene <- "Bbox1"

test <- subset(N.annot, Gene.ID == jgene)
promoterids <- as.character(unique(test$saeedid))
sitecounts <- subset(N.long, promoterid %in% promoterids)

m <- dcast(data = sitecounts, formula = promoterid ~ motif, value.var = "sitecount")
rownames(m) <- m$promoterid
m$promoterid <- NULL
m <- as.matrix(m)
m <- m[, which(colSums(m) > 0.25)]

par(mar=c(5.1, 12, 4.1, 2.1))
my_palette <- colorRampPalette(c("white", "black"))(n = 299)
heatmap.2(t(m), density.info = "density", trace = "none", margins = c(1, 14), main = jgene, col = my_palette, cexRow=0.75)


