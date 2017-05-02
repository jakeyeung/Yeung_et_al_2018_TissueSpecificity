GetSuffix <- function(jweight, use.sql, jmodstr, jcutoffstr, include.promoter=FALSE){
  if (!include.promoter){
    suffix <- paste0(".weight.", jweight, ".sql.", use.sql, ".mod.", jmodstr, ".dhscutoff.", jcutoffstr)
  } else {
    suffix <- paste0(".inclprom.", include.promoter, ".weight.", jweight, ".sql.", use.sql, ".mod.", jmodstr, ".dhscutoff.", jcutoffstr)
  }
  return(suffix)
}

GetESubDir <- function(do.center, jmodstr, jweight){
  E.subdir <- paste0("centered.", do.center, ".mod.", jmodstr, ".weightcutoff.", jweight)
}