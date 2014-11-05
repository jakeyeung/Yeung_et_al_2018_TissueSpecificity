# Jake Yeung
# Nov 5 2015
# Functions to handle sample names

ShortenSampNames <- function(long.names, show="tissue") {
  # Function to shorten sample names of a certain form. 
  # User options to determine how much to shorten
  #
  # INPUT: 
  # [GSM1321058_BS58_MoGene1.0ST.CEL, GSM1321058_BS59_MoGene1.0ST.CEL ..., ]
  # 
  # OUTPUT:
  # [BS58, BS59, ..., ]  if show == "tissue.time"
  # or
  # [BS, BS, ..., ]  if show == "tissue"
  # or
  # [58, 59, ..., ] if show == "time"
  # 
  # Whether you want the number afterwards is user defined.
  
  long.names <- strsplit(long.names, '_')
  short.names <- sapply(long.names, function(x) {
    sample.id <- x[2]    # sample name with time point
    
    # THREE WAYS OF SHOWING SAMPLES NAMES:
    
    if (show == "tissue.time") {
      # 1.
      # show sample name with timepoint e.g. Liver34
      samp.name <- sample.id 
    }
    else if (show == "tissue") {
      # 2.
      # show opnly tissue component in sample name e.g. "Adr", "Liver"
      samp.name <- substring(sample.id, 1, nchar(sample.id) - 2)      
    }
    else if (show == "time") {
      # 3.
      # show only time component in sample name e.g. 34
      samp.name <- substring(sample.id, nchar(sample.id) - 1, nchar(sample.id)) 
    }
  })
  return(short.names)
}

