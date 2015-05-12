# Load Kallisto abundance estimates and try to find first exons from there
# 2015-05-12

library(dplyr)

# Functions ---------------------------------------------------------------

source("scripts/functions/GetTissueTimes.R")
source("scripts/functions/LoadArrayRnaSeq.R")

LoadKallisto <- function(path.kallisto){
  source("scripts/functions/ConvertRNASeqTissueNamesToArray.R")
  if (missing(path.kallisto)){
    path.kallisto <- "data/alternative_exon_usage//abundance.merged.annotated.sorted.txt"
  }
  dat.kallisto <- read.table(path.kallisto, header = TRUE)
  
  # BEGIN: break matrix into a metadata and tpm estimates: makes converting to long easier
  dat.meta <- dat.kallisto[, c("chromo", "start", "end", "gene_name", "strand")]
  rownames(dat.meta) <- dat.kallisto$target_id
  
  colnames.meta <- c("chromo", "start", "end", "gene_name", "strand", "target_id")
  colnames.tpm <- colnames(dat.kallisto)[which(!colnames(dat.kallisto) %in% colnames.meta)]
  
  dat.tpm <- dat.kallisto[, colnames.tpm]
  rownames(dat.tpm) <- dat.kallisto$target_id
  # END: break matrix into a metadata and tpm estimates
  
  tissues <- sapply(colnames(dat.tpm), function(s) strsplit(s, '_')[[1]][[1]])
  tissues <- ConvertRNASeqTissueNamesToArray(tissues)
  times <- GetTimes(colnames(dat.tpm), get_unique=FALSE)
  
  tpm.long <- data.frame(transcript_id = rep(dat.kallisto$target_id, ncol(dat.tpm)),
                         chromo = rep(dat.kallisto$chromo, ncol(dat.tpm)),
                         start = rep(dat.kallisto$start, ncol(dat.tpm)), 
                         end = rep(dat.kallisto$end, ncol(dat.tpm)), 
                         gene_name = rep(dat.kallisto$gene_name, ncol(dat.tpm)),
                         strand = rep(dat.kallisto$strand, ncol(dat.tpm)),
                         tissue = rep(tissues, each = nrow(dat.tpm)),
                         time = as.numeric(rep(times, each = nrow(dat.tpm))),
                         tpm = unlist(dat.tpm))
  return(tpm.long)
}



# Load data ---------------------------------------------------------------

tpm.long <- LoadKallisto()

head(tpm.long)

jgene <- "Ddc"
test <- subset(tpm.long, gene_name == jgene)

test.afe <- test %>%
  group_by(time, tissue) %>%
  mutate(tpm_normalized = (tpm + 1) / sum(tpm + 1))
  
ggplot(test.afe, aes(y = tpm_normalized, x = transcript_id)) + geom_boxplot() + facet_wrap(~tissue) + ggtitle(jgene)


# Calculate fractional first exon usage -----------------------------------

tpm.afe <- tpm.long %>%
  group_by(gene_name, tissue, time) %>%
  mutate(tpm_normalized = tpm / sum(tpm))


