source("scripts/functions/RemoveP2Name.R")

RemoveCommasBraces <- function(m){
  return(gsub("\\{|\\}|\\,", replacement = "\\.", m))
} 

PcToMotif <- function(v, pc, top.n = 3){
  # for PMD -> MARA analysis 
  x <- v[, pc]
  x <- x[which(x != 0)]
  motif.names <- names(sort(x, decreasing = TRUE))
  return(paste(motif.names[1:top.n], collapse = "-"))
}
