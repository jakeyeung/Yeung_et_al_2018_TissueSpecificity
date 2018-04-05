# mean vs variance when you have lots of 0s
# Nov 2014
# Jake Yeung


# Init values -------------------------------------------------------------

u.o <- sapply(x1, function(x, n){
  print(x)
  X <- c(rep(x, n), rep(0, N))
  return(list(mean=mean(X), var=var(X)))
}, n)
# plot relationship between mean and variance
plot(unlist(u.o["mean", ]), unlist(u.o["var", ]), type='l', log="xy")

n.list <- seq(1:10)
N <- 1000
x1 <- seq(1:100)
for (n in n.list){
  u.o <- sapply(x1, function(x, n){
    print(x)
    X <- c(rep(x, n), rep(0, N))
    return(list(mean=mean(X), var=var(X)))
  }, n)
  # plot relationship between mean and variance
  lines(unlist(u.o["mean", ]), unlist(u.o["var", ]), type='l')
}



