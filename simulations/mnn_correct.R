args = commandArgs(trailingOnly=TRUE)
iter = as.integer(args[1])

library(Matrix)
library(batchelor)
library(NNLM)

### read in the log-pc counts
data <- readRDS(paste0("data/iter", iter, ".rds"))
norm.dat <- data$Y

### specify the patient identity information for all cells
metadata <- data.frame(group = as.factor(c(rep(1, 400), rep(2, 400), rep(3, 400), rep(4, 400), rep(5, 400), rep(6, 400), rep(7, 400), rep(8, 400))))
rownames(metadata) <- rownames(norm.dat)

### apply MNN Correct to log-pc counts
my_mnn <- mnnCorrect(t(norm.dat), batch = metadata$group)
my_mnn_integrated <- t(assay(my_mnn))

### apply NMF to the batch corrected data
my_mnn_integrated_data <- my_mnn_integrated - min(my_mnn_integrated)
fit.nmf.mnn <- nnmf(my_mnn_integrated_data, k=12, method="scd", verbose=1, rel.tol=1e-8)
saveRDS(fit.nmf.mnn, paste0("output/iter", iter, "_mnn_correct.rds"))
