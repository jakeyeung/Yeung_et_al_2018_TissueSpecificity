# Jake Yeung
# plot_clock_genes.R
# plot known clock genes and compare how they run across different tissues

# for functions: only if needed
# functions.dir <- 'scripts/functions'


# Load data ---------------------------------------------------------------


# define dirs
data_dir <- "data"
fname <- "hogenesch_2014_rma.genenames.txt"  # has duplicate gene names

# load data
data_path <- file.path(data_dir, fname)
print(paste("Reading data from,", data_path, "May take a few a minutes."))
dat <- read.table(data_path, header=TRUE)
print("Read data to memory.")

# If needed handle duplicate row names ------------------------------------

rownames(dat) <- make.names(dat$gene, unique=TRUE)

drop.cols <- c("gene")
dat <- dat[, !(names(dat) %in% drop.cols)]

Peek(dat)  # expect gene names as row names, tissues in columns

# Get Colnames ------------------------------------------------------------

colnames(dat) <- ShortenSampNames(colnames(dat), show="tissue.time")
dat.colnames <- colnames(dat)  # in case I need it later
# columns are Condition1Time1, C1T2, ... C1T24, C2T1, .... ConditionNTimeN
dat.tissuenames <- dat.colnames[seq(1, 288, 24)]  # grab every 24th column
# tissue names are Adr18... remove the 18 from tissue names
dat.tissuenames <- unname(sapply(dat.tissuenames, function(x){
  return(substr(x, 1, nchar(x) - 2))
}))

Peek(dat)  # sample names are more readable...


# Plot 24 hour rhythms ----------------------------------------------------

genes <- list("Wee1", "Per2", )
plot()
