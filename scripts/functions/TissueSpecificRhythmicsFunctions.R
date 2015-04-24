Is.TissueSpecific <- function(pval, amp, min.pval = 1e-5, max.pval = 0.05, max.amp = 0.5, min.amp = 0.1){
  # Check that min pval and max amplitude passes cutoff
  # And max pval and min amplitude passes cutoff
  
  # first two checks: if there exists a rhythmic gene
  # last two checks: if there exists a non-rhythmic gene.
  if (min(pval) < min.pval & max(amp) > max.amp & 
        max(pval) > max.pval & min(amp) < min.amp){
    return(TRUE)
  } else {
    return(FALSE)
  }
}

Is.RhythmicAcrossTissues <- function(pval, amp, cutoff.pval = 1e-5, cutoff.amp = 0.5){
  if (max(pval) < cutoff.pval & min(amp) > cutoff.amp){
    return(TRUE)
  } else {
    return(FALSE)
  }
}
