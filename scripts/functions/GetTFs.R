GetTFs <- function(){
  # Vector containing gene names (may be comma separated), get gene list
  # 
  # Args:
  # tf_vector: vector containing gene names (amy be comma separated)
  # with promoter (obtained from get_TFs_from_associations.py script)
  
  # define dirs
  data.dir <- "data"
  tf.fname <- "motifs_and_TFs.list"
  tf.path <- file.path(data.dir, tf.fname)
  
  tf.mat <- read.table(tf.path, header=FALSE, row.names = 1, sep='\t')
  
  tf_vector <- as.vector(tf.mat[, 1])
  
  genes.list <- strsplit(tf_vector, split=",")
  
  # flatten list, return uniques
  genes.list <- unique(unlist(genes.list))
  
  return(genes.list)  
}