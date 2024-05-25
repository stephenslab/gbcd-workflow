### this script shows the analyses performed on the simulated data, and reproduces the relevant figures in the paper 

### load in the required packages
library(Matrix)
library(Seurat)
library(ggplot2)
library(cowplot)
library(patchwork)
library(RColorBrewer)
library(pheatmap)
library(gridExtra)
library(fastTopics)


############################# visualize one replicate simulated dataset, as presented in the paper (Fig. 1B) #################################
### load in the simulated single cell data from this replicate
iter <- 2
dat <- readRDS(paste0("data/iter", iter, ".rds"))
X <- dat$X
L <- t(t(dat$L)/apply(dat$L, 2, max))

### create Seurat object
rownames(X) <- paste0("cell", 1:nrow(X))
colnames(X) <- paste0("gene", 1:ncol(X))
metadata <- data.frame(continuous = L[, 3], subtype = c(rep(1, 1600), rep(2, 1600)),
                       patient = c(rep(1, 400), rep(2, 400), rep(3, 400), rep(4, 400), rep(5, 400), rep(6, 400), rep(7, 400), rep(8, 400)))
rownames(metadata) <- rownames(X)

### log transform and normalize raw count data
pdac <- CreateSeuratObject(counts = t(X), project = "PDAC", meta.data = metadata)
pdac <- NormalizeData(pdac, normalization.method = "LogNormalize", scale.factor = median(rowSums(X)))

### identification of highly variable features 
pdac <- FindVariableFeatures(pdac, selection.method = "vst", nfeatures = 5000)

### scale the data
all.genes <- rownames(pdac)
pdac <- ScaleData(pdac, features = all.genes)

### perform PCA
pdac <- RunPCA(pdac, features = VariableFeatures(object = pdac), npcs = 50, verbose = FALSE)

### run tSNE
pdac <- RunTSNE(pdac, dims = 1:50)

### plot tSNE and color cells by membership of continuous GEP
p1.tsne <- FeaturePlot(pdac, features="continuous", label = FALSE, repel = TRUE) + 
  scale_color_gradient2(low = 'gray80', mid = 'blue', high = 'red', midpoint = 0.5)

### plot tSNE and color cells by subtype
p2.tsne <- DimPlot(pdac, label = FALSE, reduction = "tsne", group.by = "subtype", cols = c("lightskyblue", "magenta")) 

### plot tSNE and color cells by patient of origin
subject_col <- c(brewer.pal(n=9, name="Set1")[-9])
p3.tsne <- DimPlot(pdac, label = TRUE, reduction = "tsne", group.by = "patient", cols = subject_col) 



################## plot the GEP membership estimates of all methods for this replicate dataset (Supplementary Fig. 1 and 2) ###################
### define the color map for GEP memberships
cols <- colorRampPalette(c("gray96", "red"))(99)
brks <- seq(0, 1, length=100)

### plot the true GEP memberships
L.order <- L[, c(3,1,2,4:11)]
p.truth <- pheatmap(L.order, show_rownames = FALSE, cluster_rows = FALSE, cluster_cols = FALSE, color = cols, breaks = brks, main = "ground truth",
                    labels_col = c("continuous", paste0("subtype", 1:2), paste0("patient", 1:8)), angle_col = 90, fontsize_col = 8)

### load in the topic model fit to UMI counts
fit.tm <- readRDS(paste0("output/iter", iter, "_topic_model.rds"))
L1 <- t(t(fit.tm$L)/apply(fit.tm$L, 2, max))

### load in the NMF fit to log-pc counts
fit.nmf <- readRDS(paste0("output/iter", iter, "_nmf.rds"))
L2 <- t(t(fit.nmf$W)/apply(fit.nmf$W, 2, max))

### load in the patient-by-patient NMF fit to log-pc counts, and find highly correlated programs across patients by hierarchical clustering
res <- readRDS(paste0("output/iter", iter, "_nmf_by_patient.rds"))
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

### calculate the jaccard index between the top 50 genes for each pair of patient-specific NMF programs
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
p.cur <- pheatmap(jaccard, main="GEP signatures for programs learned by patient-by-patient nmf with k=3")

### concatenate the patient-specific NMF program memberships from the cluster of highly correlated programs across patients
programs <- sort(colnames(jaccard)[p.cur$tree_row$order][9:16])
L.sub <- L.all[, programs]
F.sub <- F.all[, programs]
L.sub <- t(t(L.sub)*apply(F.sub, 2, max))
loading <- c()
for(i in 1:8){
  loading <- c(loading, L.sub[, programs[i]])
}
L3 <- data.frame(k1=loading/max(loading))

### load in the LIGER fit to UMI counts
L4 <- readRDS(paste0("output/iter", iter, "_liger.rds"))
L4 <- t(t(L4)/apply(L4, 2, max))

### load in the CCA fit to log-pc counts implemented using Seurat v3
fit.seurat <- readRDS(paste0("output/iter", iter, "_seurat.rds"))
L5 <- t(t(fit.seurat$W)/apply(fit.seurat$W, 2, max))

### load in the MNN correct fit to log-pc counts
fit.mnn <- readRDS(paste0("output/iter", iter, "_mnn_correct.rds"))
L6 <- t(t(fit.mnn$W)/apply(fit.mnn$W, 2, max))

### load in the Conos fit to log-pc counts
fit.conos <- readRDS(paste0("output/iter", iter, "_conos.rds"))
L7 <- t(t(fit.conos$W)/apply(fit.conos$W, 2, max))

### load in the EB-SNMF fit to log-pc counts with point exponential prior on L
fit.snmf <- readRDS(paste0("output/iter", iter, "_point_exponential_snmf.rds"))
L8 <- t(t(fit.snmf$L_pm)/apply(fit.snmf$L_pm, 2, max))

### load in the EB-SNMF fit to log-pc counts with generalized binary prior on L
fit.snmf <- readRDS(paste0("output/iter", iter, "_gb_snmf.rds"))
L9 <- t(t(fit.snmf$L_pm)/apply(fit.snmf$L_pm, 2, max))

### load in the EBMF fit to log-pc counts based on covariance decomposition, with point exponential prior on L
fit.cd <- readRDS(paste0("output/iter", iter, "_point_exponential_cd.rds"))
L10 <- fit.cd$L

### load in the EBMF fit to log-pc counts based on covariance decomposition, with generalized binary prior on L
fit.gbcd <- readRDS(paste0("output/iter", iter, "_gbcd.rds"))
L11 <- fit.gbcd$L

### load in the PCA fit to log-pc counts 
fit.pca <- readRDS(paste0("output/iter", iter, "_pca_K12.rds"))
L12 <- fit.pca$u

### plot GEP membership estimates for all methods, but remove the programs that are highly correlated with cellular detection rates
loadings <- list(L1, L2, L4, L5, L6, L7, L8, L9, L10, L11)
plist <- list(NULL)
cdr <- rowMeans(X!=0)
methods <- c("NMF, combined (UMI counts)", "NMF, combined (log-pc counts)", "LIGER", "CCA (Seurat)", "MNN Correct", "Conos",
             "EB-SNMF, point-exp prior", "EB-SNMF, GB prior", "CD, point-exp prior", "GBCD", "Patient-by-patient NMF", "PCA")

for(i in 1:10){
  L.tmp <- matrix(0, nrow=nrow(L.order), ncol=0)
  k.idx <- c()
  k.cdr <- abs(cor(cdr, loadings[[i]])) > 0.7
  loadings[[i]] <- loadings[[i]][, !k.cdr]
  for(k in 1:ncol(L.order)){
    corr <- cor(L.order[,k], loadings[[i]])
    if(max(corr) > 0.6){
      k.idx <- c(k.idx, which.max(corr))
      L.tmp <- cbind(L.tmp, loadings[[i]][, which.max(corr)])
    }
  }
  if(length(k.idx) > 0)
    L.tmp <- cbind(L.tmp, loadings[[i]][, -k.idx])
  else
    L.tmp <- loadings[[i]]
  plist[[i]] <- pheatmap(L.tmp, show_rownames = FALSE, cluster_rows = FALSE, cluster_cols = FALSE, 
                         labels_col = 1:ncol(L.tmp), angle_col = 0, fontsize = 8, color = cols, breaks = brks, main = methods[i])
}

plist[[11]] <- pheatmap(L3, show_rownames = FALSE, cluster_rows = FALSE, cluster_cols = FALSE, 
                        labels_col = 1:ncol(L3), angle_col = 0, fontsize = 8, color = cols, breaks = brks, main = methods[11])

### plot all panels together
plot_grid(p.truth$gtable, plist[[1]]$gtable, plist[[2]]$gtable, plist[[3]]$gtable, plist[[4]]$gtable, plist[[5]]$gtable, 
          plist[[6]]$gtable, plist[[7]]$gtable, plist[[8]]$gtable, plist[[9]]$gtable, plist[[10]]$gtable, plist[[11]]$gtable, ncol = 3)


### plot the true versus estimated membership of the continuous GEP (Supplementary Fig. 2)
cont <- data.frame(truth=L.order[, 1], m1=L1[, which.max(cor(L1, L.order[, 1]))], m2=L2[, which.max(cor(L2, L.order[, 1]))], m3=L3[, 1],
                   m4=L4[, which.max(cor(L4, L.order[, 1]))], m5=L5[, which.max(cor(L5, L.order[, 1]))], m6=L6[, which.max(cor(L6, L.order[, 1]))], 
                   m7=L7[, which.max(cor(L7, L.order[, 1]))], m8=L8[, which.max(cor(L8, L.order[, 1]))], m9=L9[, which.max(cor(L9, L.order[, 1]))], 
                   m10=L10[, which.max(cor(L10, L.order[, 1]))], m11=L11[, which.max(cor(L11, L.order[, 1]))], 
                   m12 = L12[, which.max(cor(L12, L.order[, 1]))])
cont <- t(t(cont)/apply(cont, 2, max))
cont <- as.data.frame(cont)

p1 <- ggplot(cont, aes(x = truth, y = m1)) + geom_point(size = 1) + 
  labs(x = "True GEP membership", y = "Estimated GEP membership") + xlim(0, 1) + ylim(0, 1) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12)) +      
  annotate(geom="text", x=0.2, y=0.9, label=paste0("R=", round(cor(cont$truth, cont$m1), 2)), size=6) + ggtitle(methods[1]) 
p2 <- ggplot(cont, aes(x = truth, y = m2)) + geom_point(size = 1) + 
  labs(x = "True GEP membership", y = "Estimated GEP membership") + xlim(0, 1) + ylim(0, 1) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12)) +   
  annotate(geom="text", x=0.2, y=0.9, label=paste0("R=", round(cor(cont$truth, cont$m2), 2)), size=6) + ggtitle(methods[2]) 
p3 <- ggplot(cont, aes(x = truth, y = m3)) + geom_point(size = 1) + 
  labs(x = "True GEP membership", y = "Estimated GEP membership") + xlim(0, 1) + ylim(0, 1) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12)) +   
  annotate(geom="text", x=0.2, y=0.9, label=paste0("R=", round(cor(cont$truth, cont$m3), 2)), size=6) + ggtitle(methods[11]) 
p4 <- ggplot(cont, aes(x = truth, y = m4)) + geom_point(size = 1) + 
  labs(x = "True GEP membership", y = "Estimated GEP membership") + xlim(0, 1) + ylim(0, 1) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12)) +   
  annotate(geom="text", x=0.2, y=0.9, label=paste0("R=", round(cor(cont$truth, cont$m4), 2)), size=6) + ggtitle(methods[3]) 
p5 <- ggplot(cont, aes(x = truth, y = m5)) + geom_point(size = 1) + 
  labs(x = "True GEP membership", y = "Estimated GEP membership") + xlim(0, 1) + ylim(0, 1) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12)) +
  annotate(geom="text", x=0.2, y=0.9, label=paste0("R=", round(cor(cont$truth, cont$m5), 2)), size=6) + ggtitle(methods[4]) 
p6 <- ggplot(cont, aes(x = truth, y = m6)) + geom_point(size = 1) + 
  labs(x = "True GEP membership", y = "Estimated GEP membership") + xlim(0, 1) + ylim(0, 1) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12)) +   
  annotate(geom="text", x=0.2, y=0.9, label=paste0("R=", round(cor(cont$truth, cont$m6), 2)), size=6) + ggtitle(methods[5]) 
p7 <- ggplot(cont, aes(x = truth, y = m7)) + geom_point(size = 1) + 
  labs(x = "True GEP membership", y = "Estimated GEP membership") + xlim(0, 1) + ylim(0, 1) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12)) +   
  annotate(geom="text", x=0.2, y=0.9, label=paste0("R=", round(cor(cont$truth, cont$m7), 2)), size=6) + ggtitle(methods[6]) 
p8 <- ggplot(cont, aes(x = truth, y = m8)) + geom_point(size = 1) + 
  labs(x = "True GEP membership", y = "Estimated GEP membership") + xlim(0, 1) + ylim(0, 1) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12)) +   
  annotate(geom="text", x=0.2, y=0.9, label=paste0("R=", round(cor(cont$truth, cont$m8), 2)), size=6) + ggtitle(methods[7]) 
p9 <- ggplot(cont, aes(x = truth, y = m9)) + geom_point(size = 1) + 
  labs(x = "True GEP membership", y = "Estimated GEP membership") + xlim(0, 1) + ylim(0, 1) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12)) +   
  annotate(geom="text", x=0.2, y=0.9, label=paste0("R=", round(cor(cont$truth, cont$m9), 2)), size=6) + ggtitle(methods[8]) 
p10 <- ggplot(cont, aes(x = truth, y = m10)) + geom_point(size = 1) + 
  labs(x = "True GEP membership", y = "Estimated GEP membership") + xlim(0, 1) + ylim(0, 1) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12)) +   
  annotate(geom="text", x=0.2, y=0.9, label=paste0("R=", round(cor(cont$truth, cont$m10), 2)), size=6) + ggtitle(methods[9]) 
p11 <- ggplot(cont, aes(x = truth, y = m11)) + geom_point(size = 1) + 
  labs(x = "True GEP membership", y = "Estimated GEP membership") + xlim(0, 1) + ylim(0, 1) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12)) +   
  annotate(geom="text", x=0.2, y=0.9, label=paste0("R=", round(cor(cont$truth, cont$m11), 2)), size=6) + ggtitle(methods[10]) 
p12 <- ggplot(cont, aes(x = truth, y = m12)) + geom_point(size = 1) + 
  labs(x = "True GEP membership", y = "Estimated GEP membership") + xlim(0, 1) + ylim(0, 1) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12)) +   
  annotate(geom="text", x=0.2, y=0.9, label=paste0("R=", round(cor(cont$truth, cont$m12), 2)), size=6) + ggtitle(methods[12]) 

### plot all panels together
plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, ncol = 3) 



##################################### summarize performance of all methods across 20 replicates (Fig. 1C) ###################################
cnmf.files <- list.files("output/")

### summarize the results
res1 <- matrix(0, nrow=40, ncol=13)
res2 <- matrix(0, nrow=20, ncol=13)
res3 <- matrix(0, nrow=160, ncol=13)
idx1 <- 0
idx2 <- 0
idx3 <- 0

### iterate over replicates
for(iter in 1:20){
  
  ### load in the simulated single cell data
  dat <- readRDS(paste0("data/iter", iter, ".rds"))
  L <- dat$L
  
  #################################################### NMF methods applied to uncorrected data #######################################################
  ### load in the topic model fit to UMI counts
  fit.tm <- readRDS(paste0("output/iter", iter, "_topic_model.rds"))
  L1 <- fit.tm$L
  
  ### load in the consensus NMF fit to UMI counts
  filename <- grep(paste0("iter", iter, ".usages"), cnmf.files)
  L2 <- read.table(paste0("output/", cnmf.files[filename]))
  L2 <- L2/rowSums(L2)
  colnames(L2) <- paste0("k", 1:ncol(L2))
  
  ### load in the NMF fit to log-pc counts
  fit.nmf <- readRDS(paste0("output/iter", iter, "_nmf.rds"))
  L3 <- fit.nmf$W
  
  ### load in the patient-by-patient NMF fit to log-pc counts
  fit.nmf <- readRDS(paste0("output/iter", iter, "_nmf_by_patient_summary.rds"))
  L4 <- data.frame(k1=fit.nmf$k1)
  
  ### load in the LIGER fit
  L5 <- readRDS(paste0("output/iter", iter, "_liger.rds"))
  
  
  ############################################### NMF methods applied to to batch corrected data #################################################
  ### load in the CCA fit implemented using Seurat v3
  fit.seurat <- readRDS(paste0("output/iter", iter, "_seurat.rds"))
  L6 <- fit.seurat$W
  
  ### load in the MNN correct fit
  fit.mnn <- readRDS(paste0("output/iter", iter, "_mnn_correct.rds"))
  L7 <- fit.mnn$W
  
  ### load in the Conos fit
  fit.conos <- readRDS(paste0("output/iter", iter, "_conos.rds"))
  L8 <- fit.conos$W
  
  
  ############################################ EBMF methods applied to log-pc counts ###################################################
  ### load in the EB-SNMF fit to log-pc counts with point exponential prior on L
  fit.snmf <- readRDS(paste0("output/iter", iter, "_point_exponential_snmf.rds"))
  L9 <- fit.snmf$L_pm
  
  ### load in the EB-SNMF fit to log-pc counts with generalized binary prior on L
  fit.snmf <- readRDS(paste0("output/iter", iter, "_gb_snmf.rds"))
  L10 <- fit.snmf$L_pm
  
  ### load in the EBMF fit to log-pc counts based on covariance decomposition, with point exponential prior on L
  fit.cd <- readRDS(paste0("output/iter", iter, "_point_exponential_cd.rds"))
  L11 <- fit.cd$L
  
  ### load in the EBMF fit to log-pc counts based on covariance decomposition, with generalized binary prior on L
  fit.gbcd <- readRDS(paste0("output/iter", iter, "_gbcd.rds"))
  L12 <- fit.gbcd$L
  
  
  ############################################ PCA applied to log-pc counts ##############################################
  ### load in the PCA fit to log-pc counts
  fit.pca <- readRDS(paste0("output/iter", iter, "_pca_K12.rds"))
  L13 <- fit.pca$u
  
  
  ############################################ summarize results for all methods ##############################################
  ### estiamte the correlation between each true GEP and the estimated GEP with highest correlation among all GEPs returned by each approach
  loadings <- list(L1, L2, L3, L4, L5, L6, L7, L8, L9, L10, L11, L12, L13)

  ### iterate over subtype-related GEPs
  for(k in 1:2){
    idx1 <- idx1 + 1
    ### iterate over methods
    for(i in 1:13){
      res1[idx1, i] <- max(cor(L[,k], loadings[[i]]))
    }
  }
  
  ### iterate over the continuous GEP
  for(k in 3:3){
    idx2 <- idx2 + 1
    ### iterate over methods
    for(i in 1:13){
      res2[idx2, i] <- max(cor(L[,k], loadings[[i]]))
    }
  }
  
  ### iterate over patient-specific GEPs
  for(k in 4:11){
    idx3 <- idx3 + 1
    ### iterate over methods
    for(i in 1:13){
      res3[idx3, i] <- max(cor(L[,k], loadings[[i]]))
    }
  }
}


### summarize the estimation performance for the continuous GEP
methods <- c(rep("NMF, combined (UMI counts)", 20), rep("NMF, consensus (UMI counts)", 20), rep("NMF, combined (log-pc counts)", 20), 
             rep("Patient-by-patient NMF", 20), rep("LIGER", 20), rep("CCA (Seurat)" , 20), rep("MNN Correct", 20), rep("Conos", 20),
             rep("EB-SNMF, point-exp prior", 20), rep("EB-SNMF, GB prior", 20), rep("CD, point-exp prior", 20), rep("GBCD", 20), rep("PCA", 20))
methods <- factor(methods, 
                  levels=c("NMF, combined (UMI counts)", "NMF, consensus (UMI counts)", "NMF, combined (log-pc counts)", 
                           "Patient-by-patient NMF", "LIGER", "CCA (Seurat)", "MNN Correct", "Conos",
                           "EB-SNMF, point-exp prior", "EB-SNMF, GB prior", "CD, point-exp prior", "GBCD", "PCA"))
dat1 <- data.frame(correlation=c(res2), method=methods)

plt <- ggplot(dat1, aes(x=method, y=correlation, fill=method)) + 
  geom_boxplot(width = 0.44, size = 0.4, outlier.shape = NA, position = position_dodge(1/2)) + ylim(0, 1) +
  labs(x = "Methods", title = "Continuous GEP") + 
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=20), legend.text=element_text(size = 13.5),
        axis.text.x = element_text(size=16, hjust=1))


### summarize the estimation performance for the subtype-related GEPs
methods <- c(rep("NMF, combined (UMI counts)", 40), rep("NMF, consensus (UMI counts)", 40), rep("NMF, combined (log-pc counts)", 40), 
             rep("Patient-by-patient NMF", 40), rep("LIGER", 40), rep("CCA (Seurat)" , 40), rep("MNN Correct", 40), rep("Conos", 40),
             rep("EB-SNMF, point-exp prior", 40), rep("EB-SNMF, GB prior", 40), rep("CD, point-exp prior", 40), rep("GBCD", 40), rep("PCA", 40))
methods <- factor(methods, 
                  levels=c("NMF, combined (UMI counts)", "NMF, consensus (UMI counts)", "NMF, combined (log-pc counts)", 
                           "Patient-by-patient NMF", "LIGER", "CCA (Seurat)", "MNN Correct", "Conos",
                           "EB-SNMF, point-exp prior", "EB-SNMF, GB prior", "CD, point-exp prior", "GBCD", "PCA"))
dat2 <- data.frame(correlation=c(res1), method=methods)

plt2 <- ggplot(dat2, aes(x=method, y=correlation, fill=method)) + 
  geom_boxplot(width = 0.44, size = 0.4, outlier.shape = NA, position = position_dodge(1/2)) + ylim(0, 1) +
  labs(x = "Methods", title = "Subtype-related GEPs") + 
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=20), legend.text=element_text(size = 13.5),
        axis.text.x = element_text(size=16, hjust=1))


### summarize the estimation performance for the patient-specific GEPs
methods <- c(rep("NMF, combined (UMI counts)", 160), rep("NMF, consensus (UMI counts)", 160), rep("NMF, combined (log-pc counts)", 160), 
             rep("Patient-by-patient NMF", 160), rep("LIGER", 160), rep("CCA (Seurat)" , 160), rep("MNN Correct", 160), rep("Conos", 160),
             rep("EB-SNMF, point-exp prior", 160), rep("EB-SNMF, GB prior", 160), rep("CD, point-exp prior", 160), rep("GBCD", 160), rep("PCA", 160))
methods <- factor(methods, 
                  levels=c("NMF, combined (UMI counts)", "NMF, consensus (UMI counts)", "NMF, combined (log-pc counts)", 
                           "Patient-by-patient NMF", "LIGER", "CCA (Seurat)", "MNN Correct", "Conos",
                           "EB-SNMF, point-exp prior", "EB-SNMF, GB prior", "CD, point-exp prior", "GBCD", "PCA"))
dat3 <- data.frame(correlation=c(res3), method=methods)

plt3 <- ggplot(dat3, aes(x=method, y=correlation, fill=method)) + 
  geom_boxplot(width = 0.44, size = 0.4, outlier.shape = NA, position = position_dodge(1/2)) + ylim(0, 1) +
  labs(x = "Methods", title = "Patient-specific GEPs") + 
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=20), legend.text=element_text(size = 13.5),
        axis.text.x = element_text(size=16, hjust=1))


### plot all panels together
combined <- plt + plt2 + plt3 & theme(legend.position = "bottom")
combined + plot_layout(guides = "collect")
