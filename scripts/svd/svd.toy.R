# Toy example of SVD. Muscle genes, liver genes and "other" genes.

library(scatterplot3d)

N.liver <- 60
N.muscle <- 60
N.other <- 60
null.genes <- 500

phase.sd <- pi / 100


liver.genes <- c(complex(length.out = N.liver, modulus = rnorm(N.liver, mean = 500, sd = 10), argument = rnorm(N.liver, mean = pi / 4, sd = phase.sd)),
                 complex(length.out = N.muscle + N.other, real = rnorm((N.muscle + N.other), mean = 0, sd = 10), imaginary = rnorm((N.muscle + N.other), mean = 0, sd = 10)))
muscle.genes <- c(complex(length.out = N.liver, real = rnorm(N.liver, mean = 0, sd = 10), imaginary = rnorm(N.liver, mean = 0, sd = 10)),
                  complex(length.out = N.muscle, modulus = rnorm(N.muscle, mean = 350, sd = 10), argument = rnorm(N.muscle, mean =  3 / 4 * pi, sd = phase.sd)),
                  complex(length.out = N.other, real = rnorm(N.other, mean = 0, sd = 10), imaginary = rnorm(N.other, mean = 0, sd = 10)))
other.genes <- c(complex(length.out = N.liver + N.muscle, real = rnorm((N.liver + N.muscle), mean = 0, sd = 10), imaginary = rnorm(c(N.liver + N.muscle), mean = 0, sd = 10)),
                 complex(length.out = N.other, modulus = rnorm(N.other, mean = 25, sd = 10), argument = rnorm(N.other, mean = pi / 2, sd = phase.sd)))
both.genes <- c(complex(length.out = N.liver, modulus = rnorm(N.liver, mean = 500, sd = 10), argument = rnorm(N.liver, mean = 0, sd = phase.sd)),
                complex(length.out = N.muscle, modulus = rnorm(N.muscle, mean = 500, sd = 10), argument = rnorm(N.muscle, mean = pi, sd = phase.sd)),
                complex(length.out = N.other, modulus = rnorm(N.other, mean = 25, sd = 10), argument = rnorm(N.other, mean = pi / 2, sd = phase.sd)))

plot(liver.genes)
plot(muscle.genes)
plot(other.genes)
plot(both.genes)
                 
# liver.genes <- c(rnorm(N.liver, mean = 500, sd = 10), rnorm((N.muscle + N.other), mean = 10, sd = 10)) 
# muscle.genes <- c(rnorm(N.liver, mean = 10, sd = 10), rnorm(N.muscle, mean = 350, sd = 10), rnorm(N.other, mean = 10, sd = 10))
# other.genes <- c(rnorm((N.liver + N.muscle), mean = 0, sd = 10), rnorm(N.other, mean = 20, sd = 10))

exprs <- rbind(liver.genes, muscle.genes, other.genes, both.genes)

scatterplot3d(exprs["liver.genes", ], exprs["muscle.genes", ], exprs["other.genes", ], type = "h")

# SVD
jsvd <- svd(exprs)
jcolors <- c(rep("red", N.liver), rep("blue", N.muscle), rep("brown", N.other))

for (j.singval in c(1, 2, 3)){
  jeigen <- jsvd$v[, j.singval]
  plot(jeigen, col = jcolors, xlim = c(-max(abs(Mod(jeigen))), max(abs(Mod(jeigen)))))
  abline(v = 0)
  abline(h = 0)
}
