# plot_gene_expression_for_afe_hits.R
# Look at the top 24 AFE hits and plot their gene expression

# dat <- LoadArrayRnaSeq()

jgene <- "Gnas"
jmiso <- "14683@uc008oez.1-uc008ofa.1@uc008oer.1-uc008oes.1"
jgene <- "Map4"
jmiso <- "17758@uc009rtf.1@uc009rsy.1-uc009rsz.1-uc009rta.1-uc009rtb.1"
jgene <- "Rtn4"
jmiso <- "68585@uc007iho.1@uc007ihk.1-uc007ihl.1-uc007ihm.1"
# jgene <- "Dlgap4"
# jmiso <- "228836@uc008nnu.1@uc008nnr.1"
jgene <- "Mtus1"
jmiso <- "102103@uc009lnk.1@uc009lnm.1"

jgene <- "Acox1"
jmiso <- "69597@uc008fmg.1@uc008fmf.1-uc008fmh.1"

PlotGeneAcrossTissues(subset(dat, gene == jgene))
PlotMisoAcrossTissue(subset(summary.long, miso_id == jmiso))
