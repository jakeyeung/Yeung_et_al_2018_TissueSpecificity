source("scripts/functions/RemoveP2Name.R")

RemoveCommasBraces <- function(m){
  return(gsub("\\{|\\}|\\,", replacement = "\\.", m))
} 
