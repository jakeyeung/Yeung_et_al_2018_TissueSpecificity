# 2016-09-02
# Jake Yeung
# See how activity changes over different lambdas


rm(list=ls())


# Load --------------------------------------------------------------------

outmain.lambda <- "/home/yeung/projects/tissue-specificity/results/MARA.liver_kidney.globalLambda"
outdirs <- list.dirs(outmain.lambda, recursive = FALSE, full.names = TRUE)

outdir <- file.path(outmain.lambda, "promoters.Liver_SV129,Kidney_SV129.Liver_BmalKO,Kidney_BmalKO-Liver_SV129,Kidney_SV129.g=1001.lambda.0.0479929")
omega <- 2 * pi / 24

lapply(outdirs, function(outdir){
  indir <-  file.path(outdir, "atger_with_kidney.bugfixed")
  source("scripts/functions/LoadActivitiesLong.R")
  act.long <- LoadActivitiesLongKidneyLiver(indir, collapse.geno.tissue = TRUE, shorten.motif.name = TRUE)
  # use Liver_SV129, but should be global lambda so it shouldnt matter
  lambda.f <- file.path(indir, "Liver_SV129", "Lambda")
  lambda <- as.numeric(unlist(read.table(lambda.f)))
  act.complex <- act.long %>%
    group_by(gene, tissue) %>%
    do(ProjectToFrequency2(., omega, add.tissue=TRUE))
  act.complex <- act.complex %>%
    mutate(amp = Mod(exprs.transformed) * 4) %>%
    mutate(phase = ConvertArgToPhase(Arg(exprs.transformed), omega = omega))
  act.complex$lambda <- lambda  
  return(act.complex)
})

s.act <- SvdOnComplex(act.complex, value.var = "exprs.transformed")


max.labs <- 20
jtitle <- ""
eigens.act <- GetEigens(s.act, period = 24, comp = 1, adj.mag = TRUE, constant.amp = 4, label.n = max.labs, jtitle = jtitle, peak.to.trough = TRUE)
print(eigens.act$u.plot + ggtitle(jmod))
