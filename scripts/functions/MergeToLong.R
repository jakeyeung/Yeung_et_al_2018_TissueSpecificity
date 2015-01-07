MergeToLong <- function(normalized.array, rna.seq.exprs){
  scripts.dir <- "scripts"
  funcs.dir <- "functions"
  source(file.path(scripts.dir, funcs.dir, "GetTissueTimes.R"))
  
  tissues <- GetTissues(colnames(normalized.array))
  times.array <- GetTimes(colnames(normalized.array))
  times.rnaseq <- GetTimes(colnames(rna.seq.exprs))
  array.genes <- rownames(normalized.array)
  rnaseq.genes <- rownames(rna.seq.exprs)
  
  
  long.array <- data.frame(gene=rep(array.genes, length(tissues) * length(times.array)),
                           tissue=as.factor(rep(tissues, each=(length(array.genes) * length(times.array)))),
                           time=as.numeric(rep(times.array, each=length(array.genes))),
                           experiment=as.factor("array"),
                           exprs=unlist(normalized.array))
  
  long.rnaseq <- data.frame(gene=rep(rnaseq.genes, length(tissues) * length(times.rnaseq)),
                            tissue=rep(tissues, each=(length(rnaseq.genes) * length(times.rnaseq))),
                            time=as.numeric(rep(times.rnaseq, each=length(rnaseq.genes))),
                            experiment=as.factor("rnaseq"),
                            exprs=unlist(rna.seq.exprs))
  
  dat <- rbind(long.array, long.rnaseq)
  return(dat)
}