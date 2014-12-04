# Functions to get tissue and times: new after column fixed to array.
# Jake Yeung
# Dec 4 2014
# GetTissueTimes.R

GetTissues <- function(samp.names){
  # Samp names of form: WFAT48 (as vector)
  # return WFAT (as vector, unique)
  tissues <- unlist(lapply(samp.names, function(samp.name){
    substr(samp.name, 1, nchar(samp.name) - 2)
  }))
  return(unique(tissues))
}

GetTimes <- function(samp.names){
  # Samp names of form: WFAT48 (as vector)
  # return 48 (as vector, unique)
  times <- unlist(lapply(samp.names, function(samp.name){
    substr(samp.name, nchar(samp.name) - 1, nchar(samp.name))
  }))
  return(unique(times))
}