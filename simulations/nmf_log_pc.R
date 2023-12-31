# Set this to a number between 1 and 20 
# to load one of the simulated data sets.
iter <- 2

library(Matrix)
library(NNLM)

### read in the log-pc counts
data <- readRDS(paste0("data/iter", iter, ".rds"))
norm.dat <- data$Y

### apply NMF to the log-pc counts 
fit.nmf <- nnmf(as.matrix(norm.dat), k=12, method="scd", verbose=2, rel.tol=1e-8)
saveRDS(fit.nmf, paste0("output/iter", iter, "_nmf.rds"))
