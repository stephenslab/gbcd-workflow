# Set this to a number between 1 and 20 
# to load one of the simulated data sets.
iter <- 2

library(Matrix)
library(Seurat)
library(conos)
library(NNLM)

### read in the log-pc counts
data <- readRDS(paste0("data/iter", iter, ".rds"))
norm.dat <- data$Y

### create Seurat object
metadata <- data.frame(group = as.factor(c(rep(1, 400), rep(2, 400), rep(3, 400), rep(4, 400), rep(5, 400), rep(6, 400), rep(7, 400), rep(8, 400))))
rownames(metadata) <- rownames(norm.dat)
my_seurat <- CreateSeuratObject(counts = t(norm.dat), project = "data", meta.data = metadata)

### apply Conos to the log-pc counts
my_conos <- SplitObject(my_seurat, split.by = "group")
for (i in 1:length(my_conos)) {
  my_conos[[i]] <- FindVariableFeatures(my_conos[[i]], nfeatures = 8e3, verbose = FALSE) %>% ScaleData(verbose = FALSE) %>% RunPCA(verbose = FALSE)
}
my_conos <- Conos$new(my_conos)
my_conos$buildGraph(k = 15, k.self = 5, space = "PCA", n.odgenes = 5e3, matching.method = "mNN", 
                    metric = "angular", score.component.variance = TRUE, verbose = FALSE)
my_conos$correctGenes(n.od.genes = 8e3, verbose = FALSE)

### apply NMF to the batch corrected data
my_conos_integrated <- my_conos$expression.adj[[1]]
my_conos_integrated <- my_conos_integrated[rownames(norm.dat), ]
my_conos_integrated <- my_conos_integrated - min(my_conos_integrated)
fit.nmf.conos <- nnmf(my_conos_integrated, k=12, method="scd", verbose=1, rel.tol=1e-8)
saveRDS(fit.nmf.conos, paste0("output/iter", iter, "_conos.rds"))
