######################################### visualize HNSCC data using tSNE #########################################
### load in required packages for tSNE visualization of scRNA-seq data
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(ggrepel)
library(pheatmap)
library(gridExtra)
library(Seurat)

### load in the HNSCC data, including the library size normalized and log-transformed scRNA-seq data and annotations for cells
setwd("HNSCC")
load("hnscc.RData")

### run tSNE
hnscc <- CreateSeuratObject(counts = t(Y), project = "HNSCC", meta.data = info)
hnscc <- FindVariableFeatures(hnscc, selection.method = "vst", nfeatures = 5000)
all.genes <- rownames(hnscc)
hnscc <- ScaleData(hnscc, features = all.genes)
hnscc <- RunPCA(hnscc, features = VariableFeatures(object = hnscc), npcs = 50, verbose = FALSE)
hnscc <- RunTSNE(hnscc, dims = 1:50)

### plot tSNE colored by subject
DimPlot(hnscc, label = TRUE, repel = TRUE, pt.size = 1.5, label.size = 5, reduction = "tsne", group.by = "subject", shape.by = "subtype",
        cols = subject_col) + guides(shape = guide_legend(override.aes = list(size = 3)), ncol = 1) + theme(text = element_text(size = 20)) +
  ggtitle("") + scale_shape_manual(values = c(15, 16, 18))



##################################### apply GBCD to the HNSCC data to estimate GEP memberships and signatures #####################################
### load in required packages to run GBCD
library(Matrix)
library(ebnm)
library(flashier)
library(magrittr)
library(ashr)

### source the utility functions to fit GBCD
source("../code/fit_cov_ebnmf.R")

### fit GBCD to estimate GEP memberships and signatures
### The runtime depends on the size of the dataset being analyzed, the number of maximum GEPs and the computing environment.
### It takes about 3 hours to fit 24 GEPs for the HNSCC dataset.
fit.gbcd <- fit.cov.ebnmf(Y = Y, Kmax = 24, prior = as.ebnm.fn(prior_family = "generalized_binary", scale = 0.04), extrapolate = FALSE)
save(fit.gbcd, file = "hnscc_gbcd.RData")



######################################### examine and visualize GBCD results from the HNSCC data #########################################
### look at the GBCD results
summary(fit.gbcd)

### look at the cell membership estimates for each GEP
head(round(fit.gbcd$L, 3))

### look at the gene signatures estimates for each GEP
summary(fit.gbcd$F)
head(round(fit.gbcd$F$lfc, 3))

### add the subject and tumor subtype annotation of cells to visualization the GBCD results
anno <- data.frame(subject=info$subject, subtype=info$subtype)
rownames(anno) <- rownames(fit.gbcd$L)
anno_colors <- list(subject=subject_col, subtype=subtype_col)
cols <- colorRampPalette(c("gray96", "red"))(50)
brks <- seq(0, 1, 0.02)

### plot the annotated heatmap of GEP memberships
pheatmap(fit.gbcd$L[order(anno$subject), ], cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, 
         annotation_row = anno, annotation_colors = anno_colors, annotation_names_row = FALSE, 
         angle_col = 45, fontsize = 12, fontsize_col = 11, color = cols, breaks = brks, main = "GEP membership estimates")

### show the estimated GEP signatures and corresponding log fold change
head(round(fit.gbcd$F, 3))

### plot the volcano plot to visualize the gene signature for a given GEP
pdat <- data.frame(gene = rownames(fit.gbcd$F$lfc), 
                   lfc = fit.gbcd$F$lfc[, "GEP2"], 
                   z = abs(fit.gbcd$F$z_score[, "GEP2"]), 
                   lfsr = fit.gbcd$F$lfsr[, "GEP2"],
                   stringsAsFactors = FALSE)
pdat <- transform(pdat, lfsr = cut(lfsr, c(-1, 0.001, 0.01, 0.05, Inf)))
rows  <- with(pdat, which(!(abs(lfc) > quantile(abs(lfc), 0.996) | (z > 10))))
pdat[rows, "gene"] <- ""    
ggplot(pdat, aes(x = lfc, y = z, color = lfsr, label = gene)) + geom_point() + 
  geom_text_repel(color = "black", size = 3, segment.color = "black",
                  segment.size = 0.25, min.segment.length = 0,
                  max.overlaps = Inf, na.rm = TRUE) +
  scale_color_manual(values = c("coral", "orange", "gold", "deepskyblue")) +
  labs(x = "log-fold change", y = "|posterior z-score|") + 
  guides(colour = guide_legend(override.aes = list(size = 4))) + 
  theme(plot.title = element_text(hjust = 0.5, size = 25), axis.text = element_text(size = 20), axis.title = element_text(size = 20), 
        legend.title = element_text(size = 20), legend.text = element_text(size = 16), legend.position = "bottom") + 
  ggtitle("Volcano plot of gene signature for GEP2")