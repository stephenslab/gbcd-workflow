args = commandArgs(trailingOnly=TRUE)
iter = as.integer(args[1])

library(Matrix)
library(NNLM)

### read in the log-pc counts
data <- readRDS(paste0("data/iter", iter, ".rds"))
norm.dat <- data$Y

### apply NMF to the log-pc counts 
fit.nmf <- nnmf(as.matrix(norm.dat), k=12, method="scd", verbose=2, rel.tol=1e-8)
saveRDS(fit.nmf, paste0("output/iter", iter, "_nmf.rds"))
