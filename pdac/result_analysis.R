### this script shows the analyses performed on the PDAC data, and reproduces the relevant figures in the paper 
### note that the script for survival analysis on bulk RNA-seq data is located in the folder survival/ 

### load in the required packages
library(Matrix)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(ggrepel)
library(pheatmap)
library(gridExtra)
library(readxl)
library(dplyr)
library(Seurat)
library(stringr)
library(latex2exp)
library(openxlsx)
library(fastTopics)
library(reshape2)

### load in the sample and cohort annotation for malignant cells from all three cohorts
load("pdac_annotation.RData")

######################################### plot tsne of the single cell rna seq data (Fig. 3A) #########################################
### load in the combined UMI counts from all cohorts (we cannot share the PDAC data, because part of the data has controlled access)
load("pdac.RData")

### run tSNE
metadata <- data.frame(subject=anno$subject)
rownames(metadata) <- rownames(scdata.combined)
pdac <- CreateSeuratObject(counts = t(scdata.combined), project = "PDAC", meta.data = metadata)
pdac <- NormalizeData(pdac, normalization.method = "LogNormalize", scale.factor = median(rowSums(scdata.combined)))
pdac <- FindVariableFeatures(pdac, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(pdac)
pdac <- ScaleData(pdac, features = all.genes)
pdac <- RunPCA(pdac, features = VariableFeatures(object = pdac), npcs = 50, verbose = FALSE)
pdac <- RunTSNE(pdac, dims = 1:50)

### tSNE colored by subject
DimPlot(pdac, label = TRUE, repel = TRUE, pt.size = 1.1, label.size = 5, reduction = "tsne", group.by = "subject", cols = anno_colors$subject) + 
  guides(color = guide_legend(override.aes = list(size=2), ncol = 2)) + theme(text = element_text(size = 20)) + ggtitle("") 



##################################### apply GBCD to the PDAC data to estimate GEP memberships and signatures #####################################
### load in required packages to run GBCD
library(ebnm)
library(flashier)
library(magrittr)
library(ashr)

### load in custom functions to implement GBCD
source("../code/fit_cov_ebnmf.R")

### fit GBCD to estimate GEP memberships and signatures
### The runtime depends on the size of the dataset being analyzed, the number of maximum GEPs and the computing environment.
### It takes about 50 hours to fit 32 GEPs for the combined PDAC data containing 35,670 cells.
### Note that the "gbcd" package implements a faster version which takes about 21 hours.
fit.gbcd <- flash_fit_cov_ebnmf(Y = scdata.combined.logpc, Kmax = 32, prior = as.ebnm.fn(prior_family = "generalized_binary", scale = 0.02), 
                                maxiter = 300, extrapolate = FALSE)



######################################### plot heatmap of GEP memberships produced by gbcd (Fig. 3B) #########################################
### load in gbcd fit
fit.gbcd <- readRDS("pdac_gbcd_paper.rds")

### specify the color scheme for heatmap
cols <- colorRampPalette(c("gray96", "red"))(50)
brks <- seq(0, 1, 0.02)

### specify the indices of shared and cohort/patient-specific GEPs respectively
k.idx1 <- c(3, 30, 7, 4, 23, 12, 5, 21, 54, 15, 16, 11, 19, 28)
k.idx2 <- c(2, 56, 27, 57, 62, 58)
k.idx3 <- c(63, 40, 18, 24, 51, 32, 60, 33, 61, 25, 37, 29, 20)
k.idx <- c(k.idx1, k.idx2, k.idx3)

### determine the order of cells within each tumor by computing a 1-d embedding of GEP memberships
fit.plot <- list(L=fit.gbcd$L[, k.idx], F=fit.gbcd$L[, k.idx])
class(fit.plot) <- c("multinom_topic_model_fit")
loadings_order <- NULL
set.seed(1)
for (group in levels(anno$subject)) {
  i <- which(anno$subject == group)
  if (length(i) > 0)
    y <- drop(tsne_from_topics(select_loadings(fit.plot, i), dims = 1))
  loadings_order <- c(loadings_order, i[order(y)])
}

### plot the annotated heatmap of GEP memberships, with tumors ordered by proportion of classical-related GEP1 within each cohort
pheatmap(fit.gbcd$L[loadings_order, k.idx], cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE,
         annotation_row = anno, annotation_colors = anno_colors, annotation_names_row = FALSE, color = cols, breaks = brks,
         labels_col = 1:length(k.idx), angle_col = 0, fontsize = 7.5, fontsize_col = 11, 
         gaps_row = cumsum(table(anno$cohort))[-3], gaps_col = c(length(k.idx1), length(k.idx1)+length(k.idx2)), main = "")


### save the top driving gene lists for shared GEPs
dat.list <- list(NULL)

for(k in 1:length(k.idx1)){
  cur.genes <- rownames(fit.gbcd$F$lfc)[fit.gbcd$F$lfc[, k.idx1[k]] > pmax(quantile(fit.gbcd$F$lfc[, k.idx1[k]], 1-200/nrow(fit.gbcd$F$lfc)), log2(1.5))]
  F.tmp <- fit.gbcd$F$lfc[cur.genes,]
  cur.genes <- rownames(F.tmp)[order(F.tmp[, k.idx1[k]], decreasing = TRUE)]
  lfc <- round(fit.gbcd$F$lfc[cur.genes, k.idx1[k]], 3)
  lfsr <- fit.gbcd$F$lfsr[cur.genes, k.idx1[k]]
  zscore <- round(fit.gbcd$F$z_score[cur.genes, k.idx1[k]], 3)
  dat <- data.frame(gene=cur.genes, lfc=lfc, lfsr=lfsr, zscore=zscore)
  dat <- dat[dat$lfsr < 0.01, ]
  rownames(dat) <- 1:nrow(dat)
  dat.list[[k]] <- dat
}

names(dat.list) <- paste0("GEP", 1:length(k.idx1))
write.xlsx(dat.list, "GEP_driving_genes.xlsx")



################## determine if the GEP gene signatures show spatial structure by chromosome coordinates (Supplementary Fig. 13) #####################
library(mgcv)

### read in the gene annotation file 
anno_loc <- read.table("gencode_v19_gene_pos.txt")

### consider only the genes which are available in the gbcd fit
anno_loc <- anno_loc[anno_loc$V1 %in% rownames(fit.gbcd$F$lfc), ]

### iterate over all GEPs and all chromosomes
dev.expl.denom <- rep(0, length(k.idx))
dev.expl.numer <- rep(0, length(k.idx))
dev.expl <- matrix(NA, nrow=length(k.idx), ncol=22)
rownames(dev.expl) <- paste0("GEP", 1:length(k.idx))
colnames(dev.expl) <- paste0("chr", 1:22)

for(k in 1:length(k.idx)){
  for(i in 1:22){
    data.chr <- anno_loc[anno_loc$V2==paste0("chr", i), ]
    data.chr <- data.chr[order(data.chr$V3), ]
    data.chr$f <- fit.gbcd$F$lfc[data.chr$V1, k.idx[k]]
    data.chr$loc <- (data.chr$V3+data.chr$V4)/2
    fit.chr <- gam(f~s(loc), data=data.chr)
    dev.expl[k,i] <- summary(fit.chr)$dev.expl
    dev.expl.numer[k] <- dev.expl.numer[k] + fit.chr$deviance
    dev.expl.denom[k] <- dev.expl.denom[k] + fit.chr$null.deviance
  }
}

### calculate the proportion of deviance explained across all chromosomes
dev.expl.all <- (dev.expl.denom - dev.expl.numer)/dev.expl.denom
names(dev.expl.all) <- paste0("GEP", 1:length(k.idx))

### plot the proportion of deviance explained across all chromosomes
classes <- c(rep("shared grp1", length(k.idx1)), rep("shared grp2", length(k.idx2)-1), rep("patient-specific", length(k.idx3)))
classes <- factor(classes, levels=c("shared grp1", "shared grp2", "patient-specific"))
plt.dat <- data.frame(pve=dev.expl.all[-15], category=classes)
plt <- ggplot(plt.dat, aes(x=category, y=pve, fill=category)) + geom_boxplot(width = 0.25, size = 0.4, outlier.shape = NA) +
  scale_fill_manual(values = brewer.pal(n=12, name="Paired")[1:3]) + labs(x = "", y = TeX(r"($\rho_k$)"), title = "") + 
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 20), legend.position="none")



############################ plot heatmap of of GEP memberships produced by alternative methods (Supplementary Fig. 14-15) ############################
### pca
fit.pca <- readRDS("other/pca.rds")
L.pca <- fit.pca$u
L.pca <- t(t(L.pca)/apply(L.pca, 2, function(x){x[which.max(abs(x))]}))
rownames(L.pca) <- rownames(anno)
colnames(L.pca) <- paste0("k", 1:ncol(L.pca))
pheatmap(L.pca[loadings_order, ], cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE,
         annotation_row = anno, annotation_colors = anno_colors, annotation_names_row = FALSE, 
         color = colorRampPalette(c("blue", "gray96", "red"))(99), breaks = seq(-1, 1, length = 100),
         labels_col = 1:ncol(L.pca), angle_col = 0, fontsize = 6, fontsize_col = 11, gaps_row = cumsum(table(anno$cohort))[-3], main = "")

### consensus NMF
cnmf <- read.table("other/result.usages.k_45.dt_0_14.consensus.txt")
cnmf <- cnmf/rowSums(cnmf)
cnmf <- t(t(cnmf)/apply(cnmf, 2, max))
rownames(cnmf) <- rownames(anno)
colnames(cnmf) <- paste0("k", 1:ncol(cnmf))
pheatmap(cnmf[loadings_order, ], cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, 
         annotation_row = anno, annotation_colors = anno_colors, annotation_names_row = FALSE, color = cols, breaks = brks,
         labels_col = 1:ncol(cnmf), angle_col = 0, fontsize = 6, fontsize_col = 11, gaps_row = cumsum(table(anno$cohort))[-3], main = "")



#################### cross reference the shared GEPs identified by gbcd to the literature-derived subtype signatures (Fig. 3D) #####################
### load in signature gene sets for the subtype programs identified from Raghavan et al.
gene.sig <- read_xlsx("raghavan.xlsx", sheet=1, skip=3)
basal.r <- gene.sig$Gene[1:100]
gene.sig <- read_xlsx("raghavan.xlsx", sheet=2, skip=3)
classical.r <- gene.sig$Gene[1:100]

### load in the gene signatures for the subtype programs identified from Moffit et al.
gene.sig <- read_excel("moffit.xlsx")
basal.m <- gene.sig$symbol[order(gene.sig$F6_BasalLike, decreasing = TRUE)][1:100]
classical.m <- gene.sig$symbol[order(gene.sig$F8_Classical, decreasing = TRUE)][1:100]

### load in the gene signatures for the subtype programs identified from Chan-Seng-Yue et al.
gene.sig <- read_excel("chan_seng_yue.xlsx", skip=1, sheet=4)
classical1.c <- gene.sig$`Sig. 1 genes`[1:100]
basal1.c <- gene.sig$`Sig. 2 genes`[1:100]
classical2.c <- gene.sig$`Sig. 6 genes`[1:100]
basal2.c <- gene.sig$`Sig. 10 genes`[1:100]

### load in the gene signatures for the subtype programs identified from Hwang et al.
gene.sig <- read_excel("hwang.xlsx", skip=3, sheet=2)
classical.h <- gene.sig$`Classical-like`[1:100]
basaloid <- gene.sig$Basaloid[1:100]
squamoid <- gene.sig$Squamoid[1:100]
mesenchymal <- gene.sig$Mesenchymal[1:100]

### define the function to perform one-sided Wilcoxon rank-sum test
wilcox.test.fnc <- function(lfc, mu, sig.genes){
  ### keep only genes in the gene signature set that overlap with our analysis results
  sig.genes <- sig.genes[sig.genes %in% rownames(lfc)]
  
  ### vector of pvalues for one-sided testing
  pval <- rep(NA, ncol(lfc))
  names(pval) <- colnames(lfc)
  
  ### classify all genes into bins based on average expression levels
  expr <- rep(NA, nrow(lfc))
  names(expr) <- rownames(lfc)
  breaks <- quantile(mu, seq(0, 1, 0.02))
  for(l in 1:50){
    expr[mu >= breaks[l] & mu <= breaks[l+1]] <- l
  }
  
  ### randomly draw 100 genes from the same expression bin for each gene in the gene signature set 
  gene.idx <- c()
  for(i in 1:length(sig.genes)){
    gene.idx <- c(gene.idx, sample(which(expr==expr[sig.genes[i]]), 100))
  }
  
  ### iterate over each column of lfc
  for(k in 1:ncol(lfc)){
    pval[k] <- wilcox.test(lfc[sig.genes, k], lfc[gene.idx, k], alternative="greater")$p.value
  }
  
  ### calculate negative log pvalue 
  log.pval <- -log10(pval)
  
  return(log.pval)
}

### perform one-sided Wilcoxon rank-sum test for each signature gene set
test.basal.r <- wilcox.test.fnc(lfc=fit.gbcd$F$lfc[, k.idx1], mu=fit.gbcd$F$lfc[,1], sig.genes = basal.r)
test.classical.r <- wilcox.test.fnc(lfc=fit.gbcd$F$lfc[, k.idx1], mu=fit.gbcd$F$lfc[,1], sig.genes = classical.r)
test.basal.m <- wilcox.test.fnc(lfc=fit.gbcd$F$lfc[, k.idx1], mu=fit.gbcd$F$lfc[,1], sig.genes = basal.m)
test.classical.m <- wilcox.test.fnc(lfc=fit.gbcd$F$lfc[, k.idx1], mu=fit.gbcd$F$lfc[,1], sig.genes = classical.m)
test.classical1.c <- wilcox.test.fnc(lfc=fit.gbcd$F$lfc[, k.idx1], mu=fit.gbcd$F$lfc[,1], sig.genes = classical1.c)
test.basal1.c <- wilcox.test.fnc(lfc=fit.gbcd$F$lfc[, k.idx1], mu=fit.gbcd$F$lfc[,1], sig.genes = basal1.c)
test.classical2.c <- wilcox.test.fnc(lfc=fit.gbcd$F$lfc[, k.idx1], mu=fit.gbcd$F$lfc[,1], sig.genes = classical2.c)
test.basal2.c <- wilcox.test.fnc(lfc=fit.gbcd$F$lfc[, k.idx1], mu=fit.gbcd$F$lfc[,1], sig.genes = basal2.c)
test.classical.h <- wilcox.test.fnc(lfc=fit.gbcd$F$lfc[, k.idx1], mu=fit.gbcd$F$lfc[,1], sig.genes = classical.h)
test.basaloid <- wilcox.test.fnc(lfc=fit.gbcd$F$lfc[, k.idx1], mu=fit.gbcd$F$lfc[,1], sig.genes = basaloid)
test.squamoid <- wilcox.test.fnc(lfc=fit.gbcd$F$lfc[, k.idx1], mu=fit.gbcd$F$lfc[,1], sig.genes = squamoid)
test.mesenchymal <- wilcox.test.fnc(lfc=fit.gbcd$F$lfc[, k.idx1], mu=fit.gbcd$F$lfc[,1], sig.genes = mesenchymal)

test.sig <- cbind(test.classical.m, test.classical1.c, test.classical2.c, test.classical.r, test.classical.h,
                  test.basal1.c, test.basal.r, test.basaloid, test.mesenchymal, test.basal.m, test.basal2.c, test.squamoid)
colnames(test.sig) <- c("Classical\nMoffit", "Classical-A\nCK", "Classical-B\nCK", "Classical\nRaghavan", "Classical\nHwang",
                        "Basal-A\nCK", "Basal\nRaghavan", "Basaloid\nHwang", "Mesenchymal\nHwang", "Basal\nMoffit", "Basal-B\nCK", "Squamoid\nHwang")
rownames(test.sig) <- paste0("GEP", 1:nrow(test.sig))

### plot the dot plot for negative log pvalue
trans_fnc <- function(x){
  x <- ifelse(x < 10, 0.1, 0.1*x)
  x
}

inv_trans_fnc <- function(y){
  y <- ifelse(y < 1, 1, 10*y)
  y
}

GEP <- rep(rownames(test.sig), ncol(test.sig))
sig <- as.vector(sapply(colnames(test.sig), function(x) rep(x, nrow(test.sig))))
test.sig.plt <- data.frame(GEP=GEP, sig=sig, neglogpval=pmin(c(test.sig), 40))
test.sig.plt$GEP <- factor(test.sig.plt$GEP, levels=paste0("GEP", 1:length(k.idx1)))
test.sig.plt$sig <- factor(test.sig.plt$sig, levels=colnames(test.sig))
test.sig.plt$neglogpval2 <- test.sig.plt$neglogpval
test.sig.plt$neglogpval2[test.sig.plt$neglogpval2 > 9.5 & test.sig.plt$neglogpval2 < 10.5] <- 10
test.sig.plt$neglogpval2[test.sig.plt$neglogpval2 < 10] <- NA
p <- ggplot(test.sig.plt, aes_string(x = "GEP", y = "sig", size = "neglogpval2")) 
p <- p + geom_point() + xlab(NULL) + ylab(NULL) + theme_cowplot(font_size = 20) + guides(size = guide_legend(order = 1)) + 
  scale_size_area(name=TeX(r"($-\log_{10}$(\textit{p}-value))"), trans=scales::trans_new("trans_fnc", trans_fnc, inv_trans_fnc)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top", 
        plot.background = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank())



########################## do effect size tile plots of the top driving genes for shared GEPs identified by gbcd (Fig. 3E) ############################
### load in lfc estimates for shared GEPs identified by gbcd
F.pm <- fit.gbcd$F$lfc[, k.idx1]
colnames(F.pm) <- paste0("GEP", 1:length(k.idx1))

### make the effect size tile plots of the top 10 driving genes for each shared GEP
n <- 10
geps <- paste0("GEP", 1:length(k.idx1))
genes <- c()
for (gep in geps) {
  f <- F.pm[,gep]
  i <- head(order(f, decreasing = TRUE), n = n)
  i <- setdiff(i,genes)
  genes <- c(genes,i)
}
dat <- cbind(gene = rownames(F.pm)[genes], as.data.frame(F.pm[genes,geps]))
dat <- transform(dat, gene = factor(gene, gene))
rownames(dat) <- NULL
dat <- melt(dat)
names(dat) <- c("gene","gep","lfc")
dat <- transform(dat, gep = factor(gep,rev(geps)), lfc = lfc)
dat$lfc[dat$lfc > 5] <- 5
dat <- transform(dat, lfc = cut(lfc, breaks = c(-5, -2, -1, -0.5, 0, 0.5, 1, 2, 5)))
colors <- colorRampPalette(c("blue","white","red"))(8)
p <- ggplot(dat,aes(x = gene, y = gep, fill = lfc)) + geom_tile() + scale_fill_manual(values = colors) + 
  theme_cowplot(font_size = 8) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(x = "", y = "", fill = "LFC")



######################## do the plots characterizing in detail each shared GEP identified by gbcd (Supplementary Fig. 19) ##########################
### specify the particular GEP to plot
idx <- 14
plot.idx <- k.idx1[idx]

### plot the annotated heatmap of the GEP membership, with tumors ordered by proportion of classical-related GEP1 within each cohort
p1 <- pheatmap(fit.gbcd$L[loadings_order, plot.idx], cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE,
               annotation_row = anno, annotation_colors = anno_colors, annotation_names_row = FALSE, color = cols, breaks = brks,
               labels_col = "", fontsize = 9, gaps_row = cumsum(table(anno$cohort))[-3], main = "")

### display the gene set enrichment analysis result for the GEP signature
dat.k <- readRDS(paste0("gsea/GEP", idx, ".rds"))
dat.k <- dat.k[1:pmin(nrow(dat.k), 5),]
dat.k$cluster <- paste0("GEP", idx)
dat.k$nlog10p <- pmin(dat.k$nlog10pFishersExact, 15)
dat.k$description <- factor(dat.k$description, levels=dat.k$description)
p2 <- ggplot(dat.k, aes_string(x = "description", y = "cluster", size = "nlog10p")) + scale_x_discrete(labels = function(x) str_wrap(x, width = 20))
p2 <- p2 + geom_point(colour = "red") + xlab(NULL) + ylab(NULL) + DOSE::theme_dose(12) + scale_size_continuous(range=c(2, 6)) + 
  guides(size = guide_legend(order = 1, title = TeX(r"($-\log_{10}(P)$)"))) + ggtitle("Geneset Over-representation Analysis") +
  theme(axis.text.x = element_text(size = 15), legend.title = element_text(size = 18), plot.title = element_text(hjust = 0.5, size = 18),
        axis.text.y = element_blank(), axis.ticks.y = element_blank()) 

### make the volcano plot for the GEP signature
pdat <- data.frame(gene = rownames(fit.gbcd$F$lfc), lfc = fit.gbcd$F$lfc[, plot.idx],
                   z = abs(fit.gbcd$F$z_score[, plot.idx]), lfsr = fit.gbcd$F$lfsr[, plot.idx], stringsAsFactors = FALSE)
pdat <- transform(pdat, lfsr = cut(lfsr, c(-1, 0.001, 0.01, 0.05, Inf)))
rows  <- with(pdat, which(!(abs(lfc) > quantile(abs(lfc), 0.995) & (z > 5))))
pdat[rows, "gene"] <- ""
p3 <- ggplot(pdat, aes(x = lfc, y = z, color = lfsr, label = gene)) + geom_point() +
  geom_text_repel(color = "black", size = 3, segment.color = "black",
                  segment.size = 0.25, min.segment.length = 0,
                  max.overlaps = Inf, na.rm = TRUE) +
  scale_color_manual(values = c("coral", "orange", "gold", "deepskyblue")) +
  labs(x = "log-fold change", y = "|posterior z-score|") + 
  guides(colour = guide_legend(override.aes = list(size = 4))) + theme_cowplot(font_size = 12) + 
  theme(plot.title = element_text(hjust = 0.5, size = 25), axis.text = element_text(size = 20), axis.title = element_text(size = 20), 
        legend.title = element_text(size = 20), legend.text = element_text(size = 16), legend.position = "bottom") + ggtitle(paste0("GEP", idx))

### plot all the panels together
tmp <- cowplot::plot_grid(p3, p2, nrow = 2, rel_heights = c(2/3, 1/3))
cowplot::plot_grid(p1$gtable, tmp, ncol = 2, rel_widths = c(1/8, 7/8))
