funcs.dir <- file.path('scripts', 'functions')
source(file.path(funcs.dir, 'FourierFunctions.R'))  # for periodoigrams
source("scripts/functions/ShannonEntropy.R")

jgene <- "Defb48"
jtissues <- c("WFAT", "Liver")
jexperiments <- c("array")
for (jtissue in jtissues){
  for (jexperiment in jexperiments){
    dat.sub <- subset(dat.long, gene == jgene & tissue == jtissue & experiment == jexperiment)
    if (jexperiment == "rnaseq"){
      interval <- 6
    } else if (jexperiment == "array"){
      interval <- 2
    }
    
    freq.and.periodogram <- CalculatePeriodogram(dat.sub$exprs)
    
    freq <- freq.and.periodogram$freq[2:length(freq.and.periodogram$freq)]
    per <- freq.and.periodogram$p.scaled[2:length(freq.and.periodogram$freq)]
    
    PlotPeriodogram(freq, per, title = paste(jgene, jtissue))
    max.freqs <- FindMaxFreqs(freq, per)
    max.freqs
    max.f <- max.freqs[1]
    max.f
    max.T <- (1 / max.f) * interval
    max.T
    abline(v=max.f, col = 'blue', lwd = 2)
    text(max.f + 0.02, 0, paste0("T=", signif(max.T, digits=3), "hrs"))
    
    
    # Get Shannon entropy -----------------------------------------------------
    
    per.norm <- per / sum(per)
    
    s <- ShannonEntropy(per.norm)
    s.max <- log2(length(per.norm))
    print(jtissue)
    print(jexperiment)
    print(log(s.max / s))   
    print(per[2] / sum(per))
  }
}


