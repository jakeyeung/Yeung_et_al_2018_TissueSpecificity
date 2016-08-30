ModelStrToModel <- function(jmod){
  # Liver_SV129,Kidney_SV129.Liver_BmalKO,Kidney_BmalKO-Liver_SV129,Kidney_SV129 - >
  # c(Liver_SV129,Kidney_SV129, Liver_SV129,Kidney_SV129;Liver_BmalKO,Kidney_BmalKO)  # need some reordering magic
  jmod.long <- gsub("\\.", ";", jmod)
  jmod.long <- strsplit(jmod.long, "-")[[1]]
  # rearrange each mod so that SV129 goes before BmalKO
  jmod.long.sorted <- rep(NA, length(jmod.long))
  i <- 1
  for (j in jmod.long){
    j.sorted <- paste(sort(strsplit(j, ";")[[1]], decreasing = TRUE), collapse = ";")
    jmod.long.sorted[i] <- j.sorted
    i <- i + 1
  }
  return(jmod.long.sorted)
}
