# Set this to a number between 1 and 20 
# to load one of the simulated data sets.
iter <- 2

library(Matrix)
library(NNLM)
library(pheatmap)

### read in the log-pc counts
data <- readRDS(paste0("data/iter", iter, ".rds"))
norm.dat <- data$Y

### apply NMF with 3 programs plus an intercept term to the log-pc counts from each patient
res <- list(NULL)
for(i in 1:8){
  cur.idx <- ((i-1)*400+1):(i*400)
  cur.dat <- norm.dat[cur.idx,]
  fit.nmf <- nnmf(as.matrix(cur.dat), init=list(W0=matrix(1, nrow=nrow(cur.dat), ncol=1)), k=3, method="scd", verbose=1, rel.tol=1e-8)
  res[[i]] <- fit.nmf
}
saveRDS(res, paste0("output/iter", iter, "_nmf_by_patient.rds"))

### identify NMF programs that are shared across patients based on signature
L.all <- matrix(NA, nrow=nrow(res[[1]]$W), ncol=0)
F.all <- matrix(NA, nrow=ncol(res[[1]]$H), ncol=0)

for(i in 1:length(res)){
  L.cur <- res[[i]]$W[, c(4,1:3)]
  F.cur <- t(res[[i]]$H[c(4,1:3), ])
  colnames(L.cur) <- c(paste0("p", i, "_icpt"), paste0("p", i, "_k", 2:ncol(L.cur)))
  colnames(F.cur) <- c(paste0("p", i, "_icpt"), paste0("p", i, "_k", 2:ncol(L.cur)))
  L.all <- cbind(L.all, L.cur)
  F.all <- cbind(F.all, F.cur)
}

### calculate the jaccard index between the top 50 genes for each pair of programs
jaccard <- matrix(NA, nrow=ncol(F.all), ncol=ncol(F.all))
for(j in 1:nrow(jaccard)){
  for(k in 1:ncol(jaccard)){
    j.idx <- which(F.all[,j] > sort(F.all[,j], decreasing = TRUE)[50])
    k.idx <- which(F.all[,k] > sort(F.all[,k], decreasing = TRUE)[50])
    jaccard[j,k] <- length(intersect(j.idx, k.idx))/length(union(j.idx, k.idx))
  }
}
rownames(jaccard) <- colnames(F.all)
colnames(jaccard) <- colnames(F.all)
pheatmap(jaccard, main="GEP signatures for programs learned by patient-by-patient nmf with k=3 (plus intercept)")
