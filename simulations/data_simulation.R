library(Matrix)
library(splatter)
library(scran)
library(seqgendiff)

##################################### simulate the single cell RNA-seq data for 20 replicates ##################################################
### load in the Splatter model parameters estimated from one PDAC dataset
params <- readRDS("simparams.rds")

### define the function to normalize and log transform the UMI counts
fnc_norm <- function(X){
  ### calculate the cell-specific library size
  clusters <- quickCluster(X)
  si <- calculateSumFactors(X, clusters=clusters)
  
  ### log transform and normalize single cell count data
  norm.dat <- log(10*(median(si)*t(X)/si + 0.1))
}


### simulate single cell RNA-seq data
for(iter in 1:20){
  
  ### set the seed
  set.seed(iter)
  
  ### simulate a homoegenous population of cells
  dat <- splatSimulate(params, batchCells = 3200, seed = iter, out.prob = 0.005, lib.loc = params@lib.loc + log(2.5))
  X <- counts(dat)
  gene.info <- as.data.frame(rowData(dat))
  
  ### simulate L
  L <- matrix(0, nrow=ncol(X), ncol=11)
  L[1:1600, 1] <- 1
  L[1601:3200, 2] <- 1
  L[sample(1:nrow(L), 600, replace=FALSE), 3] <- runif(600, min=0.4, max=2)
  L[1:400, 4] <- 1
  L[401:800, 5] <- 1
  L[801:1200, 6] <- 1
  L[1201:1600, 7] <- 1
  L[1601:2000, 8] <- 1
  L[2001:2400, 9] <- 1
  L[2401:2800, 10] <- 1
  L[2801:3200, 11] <- 1
  
  ### simulate F
  F <- matrix(0, nrow=nrow(X), ncol=11)
  idx.gene <- sort(which(rowSums(X!=0) >= 300))
  F[idx.gene[1:75], 1] <- pmax(rnorm(75, log2(3), 0.5), log2(1.5))
  F[idx.gene[76:150], 2] <- pmax(rnorm(75, log2(3), 0.5), log2(1.5))
  F[idx.gene[251:500], 3] <- pmax(rnorm(250, log2(3), 0.5), log2(1.5))
  F[idx.gene[501:1000], 4] <- pmax(rnorm(500, log2(3), 0.5), log2(1.5))
  F[idx.gene[1001:1500], 5] <- pmax(rnorm(500, log2(3), 0.5), log2(1.5))
  F[idx.gene[1501:2000], 6] <- pmax(rnorm(500, log2(3), 0.5), log2(1.5))
  F[idx.gene[2001:2500], 7] <- pmax(rnorm(500, log2(3), 0.5), log2(1.5))
  F[idx.gene[2501:3000], 8] <- pmax(rnorm(500, log2(3), 0.5), log2(1.5))
  F[idx.gene[3001:3500], 9] <- pmax(rnorm(500, log2(3), 0.5), log2(1.5))
  F[idx.gene[3501:4000], 10] <- pmax(rnorm(500, log2(3), 0.5), log2(1.5))
  F[idx.gene[4001:4500], 11] <- pmax(rnorm(500, log2(3), 0.5), log2(1.5))
  F[gene.info$OutlierFactor > 1, ] <- 0
  
  ### simulate patterns of gene expression variation according to L and F using binomial thinning
  X.mod <- thin_diff(mat = as.matrix(X), design_fixed = L, coef_fixed = F)
  X.thin <- as(X.mod$mat, "sparseMatrix") 
  
  ### remove genes with very low expression levels
  idx.gene <- rowSums(X.thin!=0) >= 32
  X.thin <- X.thin[idx.gene,]
  F <- F[idx.gene,]
  colnames(X.thin) <- paste0("cell", 1:ncol(X.thin))
  rownames(X.thin) <- paste0("gene", 1:nrow(X.thin))
  rownames(L) <- colnames(X.thin)
  colnames(L) <- paste0("k", 1:ncol(L))
  rownames(F) <- rownames(X.thin)
  colnames(F) <- colnames(L)
  
  ### normalize and log transform the UMI counts
  norm.dat <- fnc_norm(X.thin)
  
  ### save the simulated data
  data <- list(X = t(X.thin), Y = norm.dat, L = L, F = F)
  saveRDS(data, file=paste0("data/iter", iter, ".rds"))
  rm(data, X, X.mod, L, F)
}
