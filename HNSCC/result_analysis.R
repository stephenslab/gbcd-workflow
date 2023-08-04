### this script shows the analyses performed on the HNSCC data, and reproduces the relevant figures in the paper 
setwd("hnscc")

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
library(fastTopics)
library(reshape2)

### load in the HNSCC data, including the library size normalized and log-transformed scRNA-seq data and annotations for cells
load("hnscc.RData")



######################################### plot tsne of the single cell rna seq data (Fig. 2A) #########################################
### run tSNE
hnscc <- CreateSeuratObject(counts = t(Y), project = "HNSCC", meta.data = info)
hnscc <- FindVariableFeatures(hnscc, selection.method = "vst", nfeatures = 5000)
all.genes <- rownames(hnscc)
hnscc <- ScaleData(hnscc, features = all.genes)
hnscc <- RunPCA(hnscc, features = VariableFeatures(object = hnscc), npcs = 50, verbose = FALSE)
hnscc <- RunTSNE(hnscc, dims = 1:50)

### plot tSNE colored by tumor
DimPlot(hnscc, label = TRUE, repel = TRUE, pt.size = 1.5, label.size = 5, reduction = "tsne", group.by = "subject", shape.by = "subtype",
        cols = subject_col) + guides(shape = guide_legend(override.aes = list(size = 3)), ncol = 1) + theme(text = element_text(size = 20)) +
  ggtitle("") + scale_shape_manual(values = c(15, 16, 18))



######################################### plot heatmap of GEP memberships produced by gbcd (Fig. 2B) #########################################
### load in gbcd fit (note this is slightly different from the gbcd fit shared in the vignette that we generated much later, 
### as we have been improving the way we implement gbcd and organize the output since we obtained the results presented in the paper)
fit.gbcd <- readRDS("hnscc_gbcd_paper.rds")

### add the corresponding cell numbers for each tumor subtype and sample
subtype <- info$subtype
subject <- info$subject
levels(subtype) <- paste0(paste(names(table(info$subtype)), table(info$subtype), sep = " ("), ")")
levels(subject) <- paste0(paste(names(table(info$subject)), table(info$subject), sep = " ("), ")")

### add the annotation of subject and subtype information for cells, and specify the color scheme for heatmap
anno <- data.frame(subject=subject, subtype=subtype)
rownames(anno) <- rownames(fit.gbcd$L)
names(subject_col) <- levels(anno$subject)
names(subtype_col) <- levels(anno$subtype)
anno_colors <- list(subject=subject_col, subtype=subtype_col)
cols <- colorRampPalette(c("gray96", "red"))(50)
brks <- seq(0, 1, 0.02)

### specify the indices of shared and patient-specific GEPs respectively
k.idx1 <- c(2, 3, 6, 17, 21, 12, 28, 47, 29, 30, 7, 13)
k.idx2 <- c(9, 20, 5, 16, 27, 8, 4)
k.idx <- c(k.idx1, k.idx2)

### determine the order of cells within each tumor by computing a 1-d embedding of GEP memberships
fit.plot <- list(L=fit.gbcd$L[, k.idx], F=fit.gbcd$L[, k.idx])
class(fit.plot) <- c("multinom_topic_model_fit")
loadings_order <- NULL
set.seed(1)
for (group in levels(subject)) {
  i <- which(subject == group)
  if (length(i) > 0)
    y <- drop(tsne_from_topics(select_loadings(fit.plot, i), dims = 1))
  loadings_order <- c(loadings_order, i[order(y)])
}

### plot the annotated heatmap of GEP memberships
pheatmap(fit.gbcd$L[loadings_order, k.idx], cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, 
         annotation_row = anno, annotation_colors = anno_colors, annotation_names_row = FALSE, color = cols, breaks = brks,
         labels_col = 1:length(k.idx), angle_col = 0, fontsize = 12, fontsize_col = 11, gaps_col = length(k.idx1), main = "")



############################# plot heatmap of of GEP memberships produced by alternative methods (Supplementary Fig. 4) #############################
### combined nmf
fit.nmf <- readRDS("other/combined_nmf.rds")
W <- t(t(fit.nmf$W)/apply(fit.nmf$W, 2, max))
W <- W[, c(24, 1:23)]
k.nmf.idx1 <- c(2, 3, 7, 8, 9, 13, 14, 16, 17, 21, 22, 24)
k.nmf.idx2 <- c(23, 12, 10, 18, 11, 6, 19, 4, 5, 20, 15)
k.nmf.idx <- c(k.nmf.idx1, k.nmf.idx2)
pheatmap(W[loadings_order, k.nmf.idx], cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, 
         annotation_row = anno, annotation_colors = anno_colors, annotation_names_row = FALSE, color = cols, breaks = brks,
         labels_col = 1:length(k.nmf.idx), angle_col = 0, fontsize = 12, fontsize_col = 11, gaps_col = length(k.nmf.idx1), main = "")

### liger
fit.liger <- readRDS("other/liger.rds")
fit.liger <- t(t(fit.liger)/apply(fit.liger, 2, max))
pheatmap(fit.liger[loadings_order, ], cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, 
         annotation_row = anno, annotation_colors = anno_colors, annotation_names_row = FALSE, color = cols, breaks = brks,
         labels_col = 1:ncol(fit.liger), angle_col = 0, fontsize = 12, fontsize_col = 11, main = "")

### seurat v3
fit.seurat <- readRDS("other/seurat.rds")
fit.seurat <- t(t(fit.seurat$W)/apply(fit.seurat$W, 2, max))
pheatmap(fit.seurat[loadings_order, ], cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, 
         annotation_row = anno, annotation_colors = anno_colors, annotation_names_row = FALSE, color = cols, breaks = brks,
         labels_col = 1:ncol(fit.seurat), angle_col = 0, fontsize = 12, fontsize_col = 11, main = "")

### mnn correct
fit.mnn <- readRDS("other/mnn_correct.rds")
fit.mnn <- t(t(fit.mnn$W)/apply(fit.mnn$W, 2, max))
pheatmap(fit.mnn[loadings_order, ], cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, 
         annotation_row = anno, annotation_colors = anno_colors, annotation_names_row = FALSE, color = cols, breaks = brks,
         labels_col = 1:ncol(fit.mnn), angle_col = 0, fontsize = 12, fontsize_col = 11, main = "")

### conos
fit.conos <- readRDS("other/conos.rds")
fit.conos <- t(t(fit.conos$W)/apply(fit.conos$W, 2, max))
pheatmap(fit.conos[loadings_order, ], cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, 
         annotation_row = anno, annotation_colors = anno_colors, annotation_names_row = FALSE, color = cols, breaks = brks,
         labels_col = 1:ncol(fit.conos), angle_col = 0, fontsize = 12, fontsize_col = 11, main = "")



###################### cross reference the shared GEPs identified by gbcd to the literature-derived meta-programs (Fig. 2B) #########################
### load in the gene signatures for the meta-programs identified by Puram et al.
gene.sig <- read_excel("puram.xlsx")
sig1 <- gene.sig$cell_cycle
sig2 <- gene.sig$`p-EMT`
sig3 <- gene.sig$epi_dif_1
sig4 <- gene.sig$epi_dif_2
sig5 <- gene.sig$stress
sig6 <- gene.sig$hypoxia

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
test.sig1 <- wilcox.test.fnc(lfc=fit.gbcd$F$lfc[, k.idx1], mu=fit.gbcd$F$lfc[,1], sig.genes = sig1)
test.sig2 <- wilcox.test.fnc(lfc=fit.gbcd$F$lfc[, k.idx1], mu=fit.gbcd$F$lfc[,1], sig.genes = sig2)
test.sig3 <- wilcox.test.fnc(lfc=fit.gbcd$F$lfc[, k.idx1], mu=fit.gbcd$F$lfc[,1], sig.genes = sig3)
test.sig4 <- wilcox.test.fnc(lfc=fit.gbcd$F$lfc[, k.idx1], mu=fit.gbcd$F$lfc[,1], sig.genes = sig4)
test.sig5 <- wilcox.test.fnc(lfc=fit.gbcd$F$lfc[, k.idx1], mu=fit.gbcd$F$lfc[,1], sig.genes = sig5)
test.sig6 <- wilcox.test.fnc(lfc=fit.gbcd$F$lfc[, k.idx1], mu=fit.gbcd$F$lfc[,1], sig.genes = sig6)
test.sig <- cbind(test.sig1, test.sig2, test.sig3, test.sig4, test.sig5, test.sig6)
colnames(test.sig) <- c("Cell\nCycle", "Partial\nEMT", "Epithelial\nDifferentiation 1", "Epithelial\nDifferentiation 2", 
                        "Cellular\nStress", "Hypoxia")
rownames(test.sig) <- paste0("GEP", 1:length(k.idx1))

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
test.sig.plt$neglogpval2[test.sig.plt$neglogpval2 > 10 & test.sig.plt$neglogpval2 < 15] <- 10
test.sig.plt$neglogpval2[test.sig.plt$neglogpval2 < 10] <- NA
p <- ggplot(test.sig.plt, aes_string(x = "GEP", y = "sig", size = "neglogpval2")) 
p <- p + geom_point() + xlab(NULL) + ylab(NULL) + theme_cowplot(font_size = 20) + guides(size = guide_legend(order = 1)) + 
  scale_size_area(name=TeX(r"($-\log_{10}$(\textit{p}-value))"), trans=scales::trans_new("trans_fnc", trans_fnc, inv_trans_fnc)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top", 
        plot.background = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank())



########################## do effect size tile plots of the top driving genes for shared GEPs identified by gbcd (Fig. 2C) ############################
### load in lfc estimates for shared GEPs identified by gbcd
F.pm <- fit.gbcd$F$lfc[, k.idx1]
colnames(F.pm) <- paste0("GEP", 1:length(k.idx1))

### make the effect size tile plots of the top 10 driving genes for each shared GEP
n <- 10
geps <- paste0("GEP", 1:length(k.idx1))
genes <- c()
for (gep in geps) {
  f <- F.pm[, gep]
  i <- head(order(f, decreasing = TRUE), n = n)
  i <- setdiff(i, genes)
  genes <- c(genes, i)
}
dat <- cbind(gene = rownames(F.pm)[genes], as.data.frame(F.pm[genes, geps]))
dat <- transform(dat, gene = factor(gene, gene))
rownames(dat) <- NULL
dat <- melt(dat)
names(dat) <- c("gene","gep","lfc")
dat <- transform(dat, gep = factor(gep, rev(geps)), lfc = lfc)
dat$lfc[dat$lfc > 10] <- 10
dat$lfc[dat$lfc < -10] <- -10
dat <- transform(dat, lfc = cut(lfc, breaks = c(-10, -5, -2, -1, -0.5, 0, 0.5, 1, 2, 5, 10)))
colors <- colorRampPalette(c("blue","white","red"))(10)
p <- ggplot(dat, aes(x = gene, y = gep, fill = lfc)) + geom_tile() + scale_fill_manual(values = colors) + 
  theme_cowplot(font_size = 8) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(x = "", y = "", fill = "LFC")



######################## do the plots characterizing in detail each shared GEP identified by gbcd (Supplementary Fig. 5) ##########################
### specify the particular GEP to plot
idx <- 1
plot.idx <- k.idx1[idx]

### plot the annotated heatmap of the GEP membership
p1 <- pheatmap(fit.gbcd$L[loadings_order, plot.idx], cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE,
               annotation_row = anno, annotation_colors = anno_colors, annotation_names_row = FALSE, color = cols, breaks = brks,
               labels_col = "", fontsize = 11, main = "")

### display the gene set enrichment analysis result for the GEP signature
dat <- readRDS(paste0("gsea/factor", plot.idx, ".rds"))
dat.k <- dat[dat$logOddsRatio > 1 & dat$nlog10pFishersExact > 5, ]
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
cowplot::plot_grid(p1$gtable, cowplot::plot_grid(p3, p2, nrow = 2, rel_heights = c(2/3, 1/3)), ncol = 2, rel_widths = c(1/6, 5/6))


