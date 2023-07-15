args = commandArgs(trailingOnly=TRUE)
iter = as.integer(args[1])

setwd("simulations")
library(Matrix)
library(ebnm)
library(flashier)
library(magrittr)

### source the utility functions to fit EB-SNMF
source("../code/fit_ebsnmf.R")

### read in the log-pc counts
data <- readRDS(paste0("data/iter", iter, ".rds"))

### fit EB-SNMF to data matrix with generalized binary prior on L and point Laplace prior on F to simultaneously estimate GEP memberships and signatures
fit.snmf <- fit.ebsnmf(Y = data$Y, Kmax = 16, prior = ebnm::ebnm_generalized_binary)
fit.snmf$sampler <- NULL
fit.snmf$flash.fit <- NULL
saveRDS(fit.snmf, paste0("output/iter", iter, "_gb_snmf.rds"))

### fit EB-SNMF to data matrix with point exponential prior on L and point Laplace prior on F to simultaneously estimate GEP memberships and signatures
fit.snmf <- fit.ebsnmf(Y = data$Y, Kmax = 16, prior = ebnm::ebnm_point_exponential)
fit.snmf$sampler <- NULL
fit.snmf$flash.fit <- NULL
saveRDS(fit.snmf, paste0("output/iter", iter, "_point_exponential_snmf.rds"))