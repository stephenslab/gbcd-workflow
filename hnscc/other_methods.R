### implement existing batch effect correction methods for single cell data, 
### including PCA, combined NMF, Harmony, Liger, Seurat v3, MNN Correct, Conos
library(Matrix)
library(Seurat)

### load in the HNSCC data, including the library size normalized and log-transformed scRNA-seq data and annotations for cells
load("hnscc.RData")

### create Seurat object
metadata <- data.frame(group = factor(info$sample.id))
rownames(metadata) <- rownames(info)
my_seurat <- CreateSeuratObject(counts = t(Y), project = "hnscc", meta.data = metadata)

### determine the number of GEPs
K <- 24


####################################### PCA ###############################################
fit.svd <- svd(Y, nu = K, nv = K)
rownames(fit.svd$u) <- rownames(Y)
rownames(fit.svd$v) <- colnames(Y)
saveRDS(fit.svd, "other/pca.rds")



####################################### combined NMF ###############################################
library(NNLM)

### run nmf with intercept
fit.nmf <- nnmf(as.matrix(Y), init=list(W0=matrix(1, nrow=nrow(Y), ncol=1)), k=K-1, method="scd", verbose=1, rel.tol=1e-8)
saveRDS(fit.nmf, "other/combined_nmf.rds")



################################################## Harmony ##################################################
library(harmony)

### apply Harmony
my_harmony <- FindVariableFeatures(my_seurat, nfeatures = 5e3, verbose = FALSE) %>% ScaleData(verbose = FALSE) %>% RunPCA(npcs = K, verbose = FALSE)
my_harmony <- RunHarmony(my_harmony, group.by.vars = "group", max.iter.harmony = 100, verbose = FALSE)
my_harmony_integrated <- Embeddings(object = my_harmony, reduction = "harmony")
sum(rownames(my_harmony_integrated)!=rownames(Y))
saveRDS(my_harmony_integrated, "other/harmony.rds")



################################################## Liger ##################################################
library(rliger)

### prepare the data list, with data from each patient as an element
my_liger <- list(NULL)
for(l in 1:nlevels(metadata$group)){
  cell.idx <- which(metadata$group==levels(metadata$group)[l])
  my_liger[[l]] <- as.matrix(t(Y[cell.idx,]))
}
names(my_liger) <- levels(metadata$group)

### apply liger 
my_liger <- createLiger(my_liger)
my_liger@norm.data <- my_liger@raw.data
my_liger <- selectGenes(my_liger)
my_liger <- scaleNotCenter(my_liger)
my_liger <- optimizeALS(my_liger, k = K)
my_liger <- quantile_norm(my_liger)
my_liger_integrated <- my_liger@H.norm[rownames(Y),]
saveRDS(my_liger_integrated, "other/liger.rds")



################################################## Seurat v3 integration ##################################################
### split the Seurat object into a list, with data from each patient as an element
my_seurat_list <- SplitObject(my_seurat, split.by = "group")
for (i in 1:length(my_seurat_list)) {
  my_seurat_list[[i]] <- FindVariableFeatures(my_seurat_list[[i]], selection.method = "vst", nfeatures = 5e3, verbose = FALSE)
}

### apply Seurat v3 integration
my_anchors <- FindIntegrationAnchors(object.list = my_seurat_list, anchor.features = 5e3, dims = 1:50, k.filter = 100, verbose = FALSE)
my_seurat_integrated <- IntegrateData(anchorset = my_anchors, dims = 1:50, k.weight = 50, verbose = FALSE)
my_seurat_integrated_data <- t(my_seurat_integrated[["integrated"]]@data[, rownames(Y)])

### apply NMF to the batch corrected data
my_seurat_integrated_data <- my_seurat_integrated_data - min(my_seurat_integrated_data)
fit.nmf.seurat <- nnmf(as.matrix(my_seurat_integrated_data), k=K, method="scd", verbose=1, rel.tol=1e-8)
saveRDS(fit.nmf.seurat, "other/seurat.rds")



################################################## MNN Correct ##################################################
library(batchelor)

### apply MNN Correct
my_mnn <- mnnCorrect(t(Y), batch = metadata$group)
my_mnn_integrated <- t(assay(my_mnn))
sum(rownames(my_mnn_integrated)!=rownames(Y))

### apply NMF to the batch corrected data
my_mnn_integrated_data <- my_mnn_integrated - min(my_mnn_integrated)
fit.nmf.mnn <- nnmf(my_mnn_integrated_data, k=K, method="scd", verbose=1, rel.tol=1e-8)
saveRDS(fit.nmf.mnn, "other/mnn_correct.rds")



################################################## Conos ##################################################
library(conos)

### apply Conos
my_conos <- SplitObject(my_seurat, split.by = "group")
for (i in 1:length(my_conos)) {
  my_conos[[i]] <- FindVariableFeatures(my_conos[[i]], nfeatures = 8e3, verbose = FALSE) %>% ScaleData(verbose = FALSE) %>% RunPCA(verbose = FALSE)
}
my_conos <- Conos$new(my_conos)
my_conos$buildGraph(k = 15, k.self = 5, space = "PCA", n.odgenes = 8e3, matching.method = "mNN", 
                    metric = "angular", score.component.variance = TRUE, verbose = FALSE)
my_conos$correctGenes(n.od.genes = 8e3, verbose = FALSE)

### apply NMF to the batch corrected data
my_conos_integrated <- my_conos$expression.adj[[1]]
my_conos_integrated <- my_conos_integrated[rownames(Y), ]
my_conos_integrated <- my_conos_integrated - min(my_conos_integrated)
fit.nmf.conos <- nnmf(my_conos_integrated, k=K, method="scd", verbose=1, rel.tol=1e-8)
saveRDS(fit.nmf.conos, "other/conos.rds")



print(sessionInfo())