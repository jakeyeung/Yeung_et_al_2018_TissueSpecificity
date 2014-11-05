# Jake Yeung
# Nov 5 2014
# test_nosvd.R
# 
# Play with rTensor library and see if decompositions match examples

library(rTensor)

vec <- c(0.9073, 0.8924, 2.1488, 
         0.7158, -0.4898, 0.3054, 
         -0.3698, 2.4288, 2.3753,
         1.7842, 1.7753, 4.2495, 
         1.6970, -1.5077, 0.3207, 
         0.0151, 4.0337, 4.7146,
         2.1236, -0.6631, 1.8260,
         -0.0740, 1.9103, 2.1335,
         1.4429, -1.7495, -0.2716)

vec <- rep(vec, 20)

A <- matrix(vec, nrow=3, ncol=180)

T <- fold(A, rs=1, cs=c(2, 3), modes=c(30, 6, 3))


# HOSVD -------------------------------------------------------------------

T.hosvd <- hosvd(T)
T.hosvd.trunc <- hosvd(T, rank=c(3,6,3))
