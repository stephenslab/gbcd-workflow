args = commandArgs(trailingOnly=TRUE)
iter = as.integer(args[1])

setwd("simulations")
library(Matrix)
library(rliger)

### read in the simulated UMI counts
dat <- readRDS(paste0("data/iter", iter, ".rds"))

### specify the patient identity information for all cells
metadata <- data.frame(group = as.factor(c(rep(1, 400), rep(2, 400), rep(3, 400), rep(4, 400), rep(5, 400), rep(6, 400), rep(7, 400), rep(8, 400))))
rownames(metadata) <- rownames(dat$X)

### prepare the data list to run liger, with data from each patient as an element
my_liger <- list(NULL)
for(l in 1:nlevels(metadata$group)){
  cell.idx <- which(metadata$group==levels(metadata$group)[l])
  my_liger[[l]] <- as.matrix(dat$X[cell.idx,])
}
names(my_liger) <- levels(metadata$group)

### apply liger 
my_liger <- createLiger(my_liger)
my_liger <- normalize(my_liger)
my_liger <- selectGenes(my_liger)
my_liger <- scaleNotCenter(my_liger)
my_liger <- optimizeALS(my_liger, k = 20)
my_liger <- quantile_norm(my_liger)
my_liger_integrated <- my_liger@H.norm
saveRDS(my_liger_integrated, paste0("output/iter", iter, "_liger.rds"))
