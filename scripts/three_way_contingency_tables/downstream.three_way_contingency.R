# 2016-10-13
# Jake Yeung

rm(list=ls())

library(dplyr)
library(ggplot2)

# Load --------------------------------------------------------------------

K <- 300
prefix <- "Robjs/three_way_cooccurence/three.way.cooccurrence.bugfixed.K."
prefix <- "Robjs/three_way_cooccurence/three.way.cooccurrence.bugfixed.nmodels.2.K."
inf <- paste0(prefix, K, ".Robj")
load(inf, v=T)

print(fits %>% arrange(pval))

# Which pairs with ONECUT?? -----------------------------------------------

jmotif <- "bHLH_family"
jmotif <- "FOXA"
jmotif <- "RORA"
jmotif <- "CUX2"
jmotif <- "ONECUT"

fits.sub <- subset(fits, grepl(jmotif, pair)) %>% arrange(pval)
print(data.frame(head(fits.sub, n = 15)))

# top guys
head(fits %>% arrange(pval))

# Plot outputs ------------------------------------------------------------

# top guys for FOXA2, ONECUT CUX2
jmotifs <- c("FOXA2", "ONECUT1.2", "CUX2")
grepstr <- paste(jmotifs, collapse = "|")

(fits.sub <- subset(fits, grepl(grepstr, pair)) %>% arrange(pval))

# contains RORA?
fits.sub$has.ROR <- sapply(fits.sub$pair, function(p) grepl("RORA", p))

fits.sub$pair <- factor(fits.sub$pair, levels = fits.sub$pair)
m <- ggplot(subset(fits.sub, pval < 0.05), aes(x = pair, y = -log10(pval), fill = has.ROR)) + geom_bar(stat = "identity") +  
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle("Top motif pairs between\nFOXA2, ONECUT, or CUX2") + 
  xlab("")
m 

# how many pairs with FOXA, ONECUT, or CUX2??
fits.sub$partner <- sapply(as.character(fits.sub$pair), function(p){
  pvec <- strsplit(p, ";")[[1]]
  partner <- pvec[!pvec %in% jmotifs]
  if (length(partner) == 0) partner <- NA
  return(partner)
})

fits.sum <- subset(fits.sub, pval < 0.05) %>%
  group_by(partner) %>%
  summarise(n.hits = length(partner)) %>%
  arrange(desc(n.hits))

fits.sum$partner <- factor(fits.sum$partner, levels = fits.sum$partner)
m.sum <- ggplot(subset(fits.sum), aes(x = partner, y = n.hits)) + geom_bar(stat = "identity")

