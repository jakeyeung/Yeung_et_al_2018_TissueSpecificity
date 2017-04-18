# 2017-04-18
# Jake Yeung

AdjustPhase <- function(phase.in, half.life, rads = FALSE, omega = 2 * pi / 24, fw.bw = "bw"){
  # Adjust phase given half.life
  gamma.mrna <- half.life / log(2)
  k.mrna <- 1 / gamma.mrna
  delay <- atan(omega / k.mrna)  # rads
  if (!rads){
    delay <- delay / omega  # rads to hours
  }
  if (fw.bw == "bw"){
    phase.adj <- phase.in - delay
  } else {
    phase.adj <- phase.in + delay
  }
  return(phase.adj)
}