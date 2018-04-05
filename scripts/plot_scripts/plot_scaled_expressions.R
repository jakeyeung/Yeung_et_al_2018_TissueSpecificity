library(dplyr)

source("scripts/functions//GetClockGenes.R")

clockgenes <- GetClockGenes()

dat.scale <- dat.long %>%
  group_by(gene, tissue, experiment) %>%
  filter(gene %in% clockgenes) %>%
  mutate(scaled = scale(exprs))

jgene <- "Arntl"
jgene2 <- "Per2"
jtitle <- jgene
jtissue <- "Mus"
jgenes <- c("Arntl", "Per2")
dat.sub <- subset(dat.scale, gene %in% jgenes & experiment == "array" & tissue == jtissue)

dat.sub$gene <- c("Bmal1", "Period2")

m <- ggplot(dat.sub, aes(x = time, y = scaled,
                     group = gene, 
                     colour = gene)) + 
  geom_line(size = 2) + 
  # facet_wrap(~tissue) +
  ggtitle("Expression of Bmal1 and Per2") + 
  ylab(label = "Normalized log expression") +
  xlab(label = "Circadian time (CT)")
m

jtissue <- c("Liver")
dat.subx <- subset(dat.scale, gene == jgene & tissue %in% jtissue)
dat.subx <- dat.subx[order(dat.subx$time), ]
dat.suby <- subset(dat.scale, gene == jgene2 & tissue %in% jtissue)
dat.suby <- dat.suby[order(dat.suby$time), ]

qplot(x = dat.suby$scaled, y = dat.subx$scaled, geom = c("point"))