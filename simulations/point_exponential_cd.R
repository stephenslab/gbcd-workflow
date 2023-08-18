# Set this to a number between 1 and 20 
# to load one of the simulated data sets.
iter <- 2

library(Matrix)
library(ebnm)
library(flashier)
library(magrittr)
library(ashr)

### source the utility functions to fit the model
source("../code/fit_cov_ebnmf.R")

### read in the log-pc counts
data <- readRDS(paste0("data/iter", iter, ".rds"))

### estimate GEP memberships and signatures by fitting EBMF based on covariance decomposition, with point exponential prior on GEP memberships
fit.cd <- flash_fit_cov_ebnmf(Y = data$Y, Kmax = 16, prior = ebnm::ebnm_point_exponential, thres = 0.9, extrapolate = FALSE)
saveRDS(fit.cd, paste0("output/iter", iter, "_point_exponential_cd.rds"))

