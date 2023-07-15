args = commandArgs(trailingOnly=TRUE)
iter = as.integer(args[1])

setwd("simulations")
library(Matrix)
library(ebnm)
library(flashier)
library(magrittr)
library(ashr)

### source the utility functions to fit the model
source("../code/fit_cov_ebnmf.R")

### read in the log-pc counts
data <- readRDS(paste0("data/iter", iter, ".rds"))

### estimate GEP memberships and signatures by fitting EBMF based on covariance decomposition, with generalized binary prior on GEP memberships
fit.gbcd <- fit.cov.ebnmf(Y = data$Y, Kmax = 16, prior = ebnm::ebnm_generalized_binary, thres = 0.9, extrapolate = FALSE)
saveRDS(fit.gbcd, paste0("output/iter", iter, "_gbcd.rds"))

