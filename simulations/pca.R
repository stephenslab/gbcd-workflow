# Set this to a number between 1 and 20 
# to load one of the simulated data sets.
iter <- 2

library(Matrix)

### read in the log-pc counts
data <- readRDS(paste0("data/iter", iter, ".rds"))
norm.dat <- data$Y

### apply PCA to the log-pc counts with 8 PCs
fit.svd <- svd(norm.dat, nu = 8, nv = 8)
rownames(fit.svd$u) <- rownames(norm.dat)
rownames(fit.svd$v) <- colnames(norm.dat)
saveRDS(fit.svd, paste0("output/iter", iter, "_pca_K8.rds"))

### apply PCA to the log-pc counts with 10 PCs
fit.svd <- svd(norm.dat, nu = 10, nv = 10)
rownames(fit.svd$u) <- rownames(norm.dat)
rownames(fit.svd$v) <- colnames(norm.dat)
saveRDS(fit.svd, paste0("output/iter", iter, "_pca_K10.rds"))

### apply PCA to the log-pc counts with 12 PCs
fit.svd <- svd(norm.dat, nu = 12, nv = 12)
rownames(fit.svd$u) <- rownames(norm.dat)
rownames(fit.svd$v) <- colnames(norm.dat)
saveRDS(fit.svd, paste0("output/iter", iter, "_pca_K12.rds"))

