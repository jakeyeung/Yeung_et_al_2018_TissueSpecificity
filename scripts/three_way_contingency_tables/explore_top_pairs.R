# 2016-11-16
# Jake Yeung

rm(list=ls())

K <- 300

load("Robjs/liver_kidney_atger_nestle/dat.long.liverkidneyWTKO.bugfixed.Robj", v=T)
load("Robjs/N.merged.pairs.Robj", v=T)

source("scripts/functions/GetTFs.R")
source("scripts/functions/HandleMotifNames.R")
source("scripts/functions/PlotUCSC.R")
source("scripts/functions/GetTFs.R")
source("scripts/functions/PlotGeneAcrossTissues.R")


# Get top pairs  ----------------------------------------------------------

# motif.pair <- c("FOXA2", "RORA")
motif.pair <- c("ONECUT1", "RORA")
motif.pair <- c("CUX2", "RORA")
motif.pair <- c("FOXA2", "RORA")
jmod <- "rhyth"  # flat or rhyth

top.bot <- c("atop", "zbottom")
N.merged$in.list <- sapply(N.merged$motif.rank, function(m) ifelse(m < K, "atop", "zbottom"))

N.merged$motif <- sapply(as.character(N.merged$motif), function(m) RemoveP2AndCommaBracesDashes(m))


# Filter for hits for both FOXA2 and RORA ---------------------------------

N.sub <- subset(N.merged, motif %in% motif.pair & model == jmod) %>%
  group_by(peak, gene) %>% 
  filter(length(in.list) == 2 & all(in.list == "atop"))

# sort peaks by ROR element


# Plot hits to output -----------------------------------------------------

# need to plot in mm9, but these coordinates are in mm10.

jsub <- dcast(N.sub, peak + gene + model ~ motif, value.var = "sitecount") %>% arrange(desc(RORA))
# jsub <- subset(N.sub, motif == "RORA") %>% arrange(desc(sitecount))
peaks <- unique(as.character(jsub$peak))

print(paste("N hits", length(peaks)))


# Calculate p-value between rhyth and flat --------------------------------

top.bot <- c("atop", "zbottom")
motif.pairs <- list(motif.pair)

# calculate atop, zbottom 2 by 2 by N tables. Where N is number of models.
N.mat.all <- dcast(subset(N.merged, motif %in% motif.pair), formula = gene + peak + model ~ motif, value.var = "in.list", fill = top.bot[2])

start <- Sys.time()
N.mat.freqs <- lapply(motif.pairs, function(motif.pair){
  motif.pair.str <- paste(motif.pair, collapse = ";")
  N.mat.freq <- N.mat.all %>%
    group_by_(.dots = c("model", motif.pair)) %>%
    summarise(freq = length(gene)) %>%
    mutate(pair = motif.pair.str)
  return(N.mat.freq)
})
print(Sys.time() - start)

N.mat.freqs <- rbindlist(N.mat.freqs)  # may not be accurate because you filter for pairs in your analysis 
colnames(N.mat.freqs) <- c("model", "motif1", "motif2", "freq", "pair")


# Write bed  --------------------------------------------------------------


bed <- data.frame(chromo = sapply(peaks, function(p) strsplit(p, ":")[[1]][[1]]),
                  startend = sapply(peaks, function(p) strsplit(p, ":")[[1]][[2]]), stringsAsFactors = FALSE)
bed$start <- sapply(bed$startend, function(s) strsplit(s, "-")[[1]][[1]], simplify = TRUE)
bed$end <- sapply(bed$startend, function(s) strsplit(s, "-")[[1]][[2]], simplify = TRUE)


pairname <- paste(motif.pair, collapse = "_")
fname.mm9 <- paste0("pairs.", pairname, ".mm9.K", K, ".bed")
fname.mm10 <- paste0("pairs.", pairname, ".mm10.K", K, ".bed")
fname.pdf <- paste0("pairs.", pairname, ".mm9.K", K, ".pdf")

write.table(subset(bed, select = c(chromo, start, end)), file = paste0("/home/shared/tissue_specificity/tissue_clock_motif_pairs/", fname.mm10), 
            quote = FALSE, append = FALSE, sep = "\t", row.names = FALSE)


# Map mm10 to mm9 ---------------------------------------------------------


cmd <- paste0("bash /home/yeung/projects/ucsc_utils/mm10_to_mm9_bed.sh /home/shared/tissue_specificity/tissue_clock_motif_pairs/", fname.mm10, 
              " /home/shared/tissue_specificity/tissue_clock_motif_pairs/unmapped.bed /home/shared/tissue_specificity/tissue_clock_motif_pairs/", fname.mm9)
system(cmd)


# Print screenshots of mm9  -----------------------------------------------

mm9f <- paste0("/home/shared/tissue_specificity/tissue_clock_motif_pairs/", fname.mm9)

(bed.mm9 <- read.table(mm9f, header = TRUE))


bedToUCSC(toPlot = bed.mm9, outpdf = paste0("plots/ucsc_motif_screenshots/", fname.pdf), 
          leftwindow = 0, rightwindow = 0, 
          theURL = "http://genome-euro.ucsc.edu/cgi-bin/hgTracks?hgsid=218772886_HoC7q8PAUF1alvtA9a4LXP0OOAcL&hgt.psOutput=on")
