# Jake Yeung
# cbind_activities.R
# After running MARA for each tissue, we get a file of activities for each tissue.
# Combine each activity from each tissue into a long data frame (makes plotting easy).
# 2015-30-03



# Get dir info ------------------------------------------------------------


act.dir <- "/home/yeung/projects/tissue-specificity/results/outputs_all_genes_MARA.swissregulon.corrected/MARA.33motifs.corrected"
out.fpath <- "/home/yeung/projects/tissue-specificity/results/outputs_all_genes_MARA.swissregulon.corrected/MARA.33motifs.corrected/activities.all"
tissues <- dir(act.dir)
times.vec <- seq(18, 64, 2)


# Loop through tissue directories -----------------------------------------

for (i in 1:length(tissues)){
  tissue <- tissues[i]
  act.path <- file.path(act.dir, tissue, "Activities")
  act <- read.table(act.path, header = FALSE, row.names = 1)
  cnames <- paste0(tissue, times.vec)
  colnames(act) <- cnames
  
  # cbind each activities dataframe to each other.
  # if first dataframe, then initialize.
  if (i == 1){
    act.all <- act
  } else {
    # check rownames match with act.all so we don't cbind stupid things
    if (all(rownames(act) == rownames(act.all))) {
      act.all  <- cbind(act.all, act)
    } else {
      warning(paste(tissue, "does not have same rownames. Skipping..."))
    }
  }
}


# Write to file -----------------------------------------------------------

write.table(act.all, file = out.fpath, quote = FALSE, sep = "\t", col.names = NA, row.names = TRUE)
