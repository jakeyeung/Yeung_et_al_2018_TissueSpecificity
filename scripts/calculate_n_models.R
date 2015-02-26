# calculate_n_models.R
# Calculate number of models from Cedric's code.

library(combinat)

for (N in seq(10)){
  n.rhythmic <- 0
  n.rhythmic.osp <- 0
  
  for (i in seq(0, N)){
    n.rhyth.temp <- nCm(N, i)
    n.rhythmic <- n.rhythmic + n.rhyth.temp
    for (k in seq(2, max(i, 2))){
      print(nCm(i, k))
      n.rhythmic.osp <- n.rhythmic.osp + nCm(i, k) * n.rhyth.temp
    }
  }
  print(paste("N:", N))
  print(paste('n.rhythmic:', n.rhythmic))
  print(paste('n.rhythmic.osp:', n.rhythmic.osp))
  print(paste('total', n.rhythmic + n.rhythmic.osp))  
}

