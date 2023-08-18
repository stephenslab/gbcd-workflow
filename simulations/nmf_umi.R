# Set this to a number between 1 and 20 
# to load one of the simulated data sets.
iter <- 2

library(Matrix)
library(fastTopics)

### read in the simulated UMI counts
data <- readRDS(paste0("data/iter", iter, ".rds"))

### fit topic model to UMI counts
fit.tm <- fit_topic_model(data$X, k = 11, numiter.main = 100, numiter.refine = 500, verbose = "progressbar")
saveRDS(fit.tm, paste0("output/iter", iter, "_topic_model.rds"))  
