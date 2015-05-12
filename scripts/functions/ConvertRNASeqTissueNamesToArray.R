ConvertRNASeqTissueNamesToArray <- function(jnames){
  array.names <- c("Adr", "Aorta", "BFAT", "BS", "Cere", "Heart", "Hypo", "Kidney", "Liver", "Lung", "Mus", "WFAT")
  rnaseq.names <- c('Adr','Aor','BFat','Bstm','Cer','Hrt','Hyp','Kid','Liv','Lun','Mus','WFat')
  names.dic <- setNames(object = array.names, nm = rnaseq.names)
  return(names.dic[jnames])
}