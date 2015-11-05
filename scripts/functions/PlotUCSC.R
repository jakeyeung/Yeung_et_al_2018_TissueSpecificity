
# Functions ---------------------------------------------------------------

####### lib ############

# Here is an R script wrote by Aaron Statham which saves UCSC to pdfs -
# you can choose which genome and tracks to display by altering the 'url' parameter. 'trackfile' is the url of a file describing the custom tracks (beds/bigwigs) to display
mergePDF <- function(output="merged.pdf", sourcefiles=c("source1.pdf","source2.pdf","source3.pdf"))
{
  # create the command string and call the command using system()
  # merging command from http://hints.macworld.com/article.php?story=2003083122212228
  command=paste("gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite",paste("-sOutputFile",output, sep="="), paste(sourcefiles, collapse=" "),sep=" ")
  try(system(command))
}


#  Reference: http://www.biostars.org/post/show/6132/batch-viewing-of-ucsc-browser-graphic/
screenshotUCSC <- function(url, trackfile, chr, start, end, filename) {
  oldpen <- options("scipen")
  options(scipen=100)
  temp <- readLines(paste(url, "&hgt.customText=", trackfile, "&position=",chr,":",start,"-",end, sep=""))
  # cat(temp,"\n")
  pdfurl <- paste("http://www.genome.ucsc.edu/trash",gsub(".*trash","",gsub(".pdf.*","",temp[grep(".pdf", temp, fixed=TRUE)][1])), ".pdf", sep="")
  # cat(pdfurl,"\n");
  options(scipen=oldpen)
  download.file(pdfurl, filename, mode="wb", quiet=TRUE)
}

bedToUCSC <- function(toPlot, outpdf, leftwindow = 2000, rightwindow = 1999, jtempdir, theURL, chromo.name="chromo", start.name="start", end.name="end"){
  # from http://onetipperday.blogspot.ch/2012/07/get-ucsc-images-for-list-of-regions-in.html
  
  # example of controling individual track
  # theURL: grab the url from the View -> PDF page of your genome browser session 
  # theURL="http://www.genome.ucsc.edu/cgi-bin/hgTracks?hgsid=451672941_fba5iOL1AK9ZwLVNqWeQrZj5fyTq&hgt.psOutput=on"
  # read regions
  # toPlot=read.table("/home/yeung/projects/tissue-specificity/data/alternative_exon_usage/abundance.merged.annotated.sorted.pythonmerged.bed", header=T)
  ## paralle version
  # library(multicore)
  # mclapply(1:nrow(toPlot), function(i) screenshotUCSC(theURL, "", as.character(toPlot$chr[i]), toPlot$start[i]-2000, toPlot$end[i]+1999, paste("region_", i, "_", toPlot$name[i],".pdf", sep="")), mc.cores=10)
  
  # anti-robot version 
  # UCSC Policy: Program-driven use of this software is limited to a maximum of one hit every 15 seconds and no more than 5,000 hits per day.
  
  if(missing(jtempdir)){
    jtempdir <- "/home/yeung/projects/tissue-specificity/plots/temp"
  }
  if(missing(theURL)){
    theURL <- "http://www.genome.ucsc.edu/cgi-bin/hgTracks?hgsid=451672941_fba5iOL1AK9ZwLVNqWeQrZj5fyTq&hgt.psOutput=on"
  }
  for(i in 1:nrow(toPlot)){
    tempname <- paste("region_", i, "_", toPlot$name[i],".pdf", sep="")
    screenshotUCSC(theURL, "", as.character(toPlot[[chromo.name]][i]), toPlot[[start.name]][i]-leftwindow, toPlot[[end.name]][i]+rightwindow, file.path(jtempdir, tempname))
    Sys.sleep(5) 
  }
  
  # merge script
  tempfiles.name <- list.files(jtempdir, pattern="region_.*pdf")
  tempfiles.path <- file.path(jtempdir, tempfiles.name)
  mergePDF(outpdf, tempfiles.path)
  rmstr <- paste0("rm ", jtempdir, "/region_.*pdf")
  try(system(rmstr))
  
}



# Main --------------------------------------------------------------------



