# svd.on.rhythmic.genes.analysis.R
# continuation from svd.on.rhythmic.genes.R
# Feb 23 2015

GetFractionOfSignal <- function(component, original){
  proj.vec <- ProjectOnVector(component, original)
  return(proj.vec / Mod(original))
}

ComposeMatrixFromComponent <- function(s, component){
  # after SVD decomposition, reconstruct a matrix from a component
  s.mat <- s$d[component] * OuterComplex(s$u[, component, drop = FALSE], t(s$v[, component, drop = FALSE]))
  return(s.mat)
}

# # Analyze contribution of component to signal -----------------------------

# Take top 100 genes from a component, see whether it contributes significantly
# to the original expression

gene <- "Myh7"

orig <- Mod(as.matrix(dat.wide[gene, ]))
proj <- matrix(data = 0, nrow = 1, ncol = 11)
for (component in 1:11){
  egene <- s$v[, component, drop = FALSE]
  esample <- s$u[, component, drop = FALSE]
  evalue <- s$d[component]
  component.mat <- evalue * OuterComplex(esample, t(egene))

  proj.vec <- ProjectOnVector(component.mat[gene, ], dat.wide[gene, ])
  print(paste('Component', component))
  print(proj.vec / orig)
  
  proj <- proj + proj.vec / orig
}

print(proj)


# Toy SVD examples --------------------------------------------------------

# Test out penalized regression
# install.packages("PMA")  # requires biocLite("impute")

max.components <- 4

# add two components together
dat.toy <- matrix(data = 0, nrow = 20973, ncol = 11)
for (component in 1:max.components){
  dat.toy <- dat.toy + ComposeMatrixFromComponent(s, component)
}

# remove zeros
dat.toy <- dat.toy[rowSums(Mod(dat.toy))>0, ]

# # Add fake.gene to dat.toy 
# fake.gene <- dat.wide["Apobec2", ]
# rownames(fake.gene) <- "fake.gene"
# dat.toy <- rbind(dat.toy, fake.gene)
# 
# # Add holes to dat.toy
holes <- list(
#   "Arntl" = "Adr",
#   "Dbp" = "Aorta",
#   "Dtx4" = "BFAT",
#   "Nr1d1" = "BS",
#   "Npas2" = "Cere",
#   "Nr1d2" = "Heart",
#   "Per2" = "Hypo",
#   "Per3" = "Kidney",
#   "Lonrf3" = "Liver",
#   "Nfil3" = "Lung",
#   "Tef" = "Mus",
#   "Dtx4" = "Adr"
)

# ADD MISSING HOLES to see effect on svd
for (gene in names(holes)){
  missing.tiss <- holes[[gene]]
  dat.toy[gene, missing.tiss] <- 0.0001  # my hole, non zero to prevent infinities
}

# svd my matrix
# s.toy <- svd(dat.toy)

library(devtools)
dev_mode(on=T)
install_github("jakeyeung/PMA")
library(PMA)
s.toy.cv <- PMD.cv(dat.toy, type = "standard", sumabss = seq(1 / sqrt(ncol(dat.toy)), 1, len = 10))
s.toy <- PMD(dat.toy, sumabsu = sqrt(nrow(dat.toy)), sumabsv = sqrt(ncol(dat.toy)), K = 7, center = FALSE)
s.toy <- PMD(dat.toy, sumabs = 1, K = 7, center = FALSE)
rownames(s.toy$u) <- rownames(dat.toy)
rownames(s.toy$v) <- colnames(dat.toy)
# when finished do:
dev_mode(on=F) #and you are back to having stable ggplot2




# plot my eigenvalues
plot(s.toy$d[1:length(s.toy$d)])

# my orig which is a superposition of components
orig <- t(Mod(as.matrix(dat.toy[gene, ])))

missing.tissues <- unname(sapply(holes, function(h) return(h[[1]])))
for (gene in names(holes)){
  for (component in 1:max.components){
    s.toy.mat <- s.toy$d[component] * OuterComplex(s.toy$u[, component, drop = FALSE], t(s.toy$v[, component, drop = FALSE]))
    frac.signal <- GetFractionOfSignal(s.toy.mat[gene, missing.tissues], dat.toy[gene, missing.tissues])
    print(paste("Component:", component, "gene", gene))
    print(frac.signal) 
  }
}

for (component in 1:max.components){
  s.toy.mat <- s.toy$d[component] * OuterComplex(s.toy$u[, component, drop = FALSE], t(s.toy$v[, component, drop = FALSE]))
  PlotComplex(s.toy$d[component] * s.toy$v[, component], labels = rownames(s.toy$v), main = paste("Component", component), add.text.plot = FALSE)

  eigensample <- s.toy$u[, component]
  jmax <- max(Mod(eigensample))
  # rotate opposite way to make it all kosher
  PlotComplex(eigensample,
              axis.min = -jmax,
              axis.max = jmax,
              labels = names(eigensample),
              col = "HSV",
              main = paste("Component:", component),
              add.text.plot = FALSE,
              jpch = 1,
              threshold = 0.5 * jmax,
              verbose = FALSE,
              rotate = jmax.arg)
  
}

