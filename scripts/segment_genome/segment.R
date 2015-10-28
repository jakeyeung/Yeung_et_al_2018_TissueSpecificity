# Segment: find genes that correlate up and down in expression by genomic proximity
# 2015-10-28


# Source ------------------------------------------------------------------

source("scripts/functions/Segment.R")



n=100
M1=make.block(rnorm(n), rnorm(12))
M2=make.block(rnorm(n), rnorm(5))
M3=make.block(rnorm(n), rnorm(5))
M4=make.block(rnorm(n), rnorm(12))
M5=make.block(rnorm(n), rnorm(5))

# M=cbind(M1,M2)
M=cbind(M1, M2, M3, M4, M5)
M=M + rnorm(length(M))/10  # noise level should be 0.01

# res=seg(M)
# print(res)

# if BIC
res=seg.dp(M, sigma2=0.1, BIC=T)
print(res)
