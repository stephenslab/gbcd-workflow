args = commandArgs(trailingOnly=TRUE)
iter = as.integer(args[1])

library(Matrix)
library(Seurat)
library(NNLM)

### read in the log-pc counts
data <- readRDS(paste0("data/iter", iter, ".rds"))
norm.dat <- data$Y

### create Seurat object
metadata <- data.frame(group = as.factor(c(rep(1, 400), rep(2, 400), rep(3, 400), rep(4, 400), rep(5, 400), rep(6, 400), rep(7, 400), rep(8, 400))))
rownames(metadata) <- rownames(norm.dat)
my_seurat <- CreateSeuratObject(counts = t(norm.dat), project = "data", meta.data = metadata)

### split the Seurat object into a list, with data from each patient as an element
my_seurat_list <- SplitObject(my_seurat, split.by = "group")
for (i in 1:length(my_seurat_list)) {
  my_seurat_list[[i]] <- FindVariableFeatures(my_seurat_list[[i]], selection.method = "vst", nfeatures = 5000, verbose = FALSE)
}

### apply Seurat v3 integration based on CCA to the log-pc counts
my_anchors <- FindIntegrationAnchors(object.list = my_seurat_list, anchor.features = 5000, dims = 1:50, verbose = FALSE)
my_seurat_integrated <- IntegrateData(anchorset = my_anchors, dims = 1:50, verbose = FALSE)
my_seurat_integrated_data <- t(my_seurat_integrated[["integrated"]]@data)

### apply NMF to the batch corrected data
my_seurat_integrated_data <- my_seurat_integrated_data - min(my_seurat_integrated_data)
fit.nmf.seurat <- nnmf(as.matrix(my_seurat_integrated_data), k=12, method="scd", verbose=1, rel.tol=1e-8)
saveRDS(fit.nmf.seurat, paste0("output/iter", iter, "_seurat.rds"))
