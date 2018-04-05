# cluster_models.R
# Jake Yeung
# 11 March 2015


# Functions ---------------------------------------------------------------

LoadAndFilterDat <- function(){
  load("~/projects/select-rhythmic-models/data/tissue_array_dat.Robj")  # dat
  conds <- c('Adr', 'Aorta', 'BFAT', 'Kidney', 'Liver', 'Lung', 'Mus', 'Heart')  # order matters
  co <- length(conds)
  t <- rep(seq(18, 64, 2), co)
  
  # filter by conditions
  conds.grep <- paste0(conds, collapse = "|")
  dat.filt <- dat[, grepl(conds.grep, colnames(dat))]
  return(dat.filt)
}


# Main --------------------------------------------------------------------

setwd('/home/yeung/projects/tissue-specificity')
dat.dir <- 'nconds_outputs/nconds_8_conds_rerun_less_cores'
load(file.path(dat.dir, 'nconds_8_tissues.matrix_output.Robj'))  # my_mat
load(file.path(dat.dir, 'nconds_8_tissues.fit_output.Robj'))  # fit
dat <- LoadAndFilterDat()

library("devtools")
dev_mode()
install("~/projects/ncond")  # use master branch
library(nconds)
library(parallel)

t <- GetTVector(dat)
n.co <- GetNConds(dat)
conds <- GetConds(dat)
period <- 24

start <- Sys.time()
fit.all_bic = parallel::mclapply(split(dat, rownames(dat)),
                                 do_all,
                                 t=t,
                                 n.co=n.co,
                                 period=period,
                                 my_mat=my_mat,
                                 conds=conds,
                                 mc.cores=10,
                                 return_min_bic=FALSE,
				 do_gc=TRUE)
print(Sys.time() - start)
save(fit.all_bic, file = 'nconds_outputs/nconds_8_conds.fit_all_bic.Robj')
# dat.fit <- InsertFitToMat(fit, dat, 8)



dev_mode(F)
