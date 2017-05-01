GetSuffix <- function(jweight, use.sql, jmodstr, jcutoffstr){
  suffix <- paste0(".weight.", jweight, ".sql.", use.sql, ".mod.", jmodstr, ".dhscutoff.", jcutoffstr)
}

GetESubDir <- function(do.center, jmodstr, jweight){
  E.subdir <- paste0("centered.", do.center, ".mod.", jmodstr, ".weightcutoff.", jweight)
}