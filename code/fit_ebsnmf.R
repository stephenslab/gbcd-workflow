### fit Empirical Bayes Semi Nonnegative Matrix Factorization (EB-SNMF) to single cell RNA-seq data, with a nonnegative prior on GEP memberships
### Y: a cell by gene matrix of normalized and log-transformed gene expression data
### Kmax: a positive integer (at least 2) specifying an upper bound for the number of GEPs to be identified
### prior: a nonnegative prior for GEP memberships, such as generalized binary prior, which must be a function defined in the ebnm package
### extrapolate: a logical indicator specifying whether to use extrapolation to accelerate backfitting GEP memberships; see the flashier package for details
### maxiter: a positive integer specifying the maximum number of iterations to backfit GEP memberships
### verbose: an integer specifying whether and how to display progress updates, as described in the flashier package
flash_fit_ebsnmf <- function(Y, Kmax, prior = ebnm::ebnm_generalized_binary, extrapolate = TRUE, maxiter = 500, verbose = 1){

  ### fit EBMF with point laplace prior to gene expression data matrix
  fit.mf <- flash_init(Y, S = 1/sqrt(nrow(Y)), var_type = 2) %>%
    flash_greedy(Kmax = Kmax, ebnm_fn = ebnm_point_laplace) %>%
    flash_backfit(maxiter = maxiter, verbose = verbose)

  ### initialize EB-SNMF fit from the EBMF fit with point laplace prior
  snmf.init <- init.ebsnmf(fit.mf, 2e1)

  ### fit EB-SNMF with a nonnegative prior on L to gene expression data matrix
  fit.snmf <- flash_init(Y, S = 1/sqrt(nrow(Y)), var_type = 2) %>%
    flash_factors_init(
      init = snmf.init,
      ebnm_fn = c(prior, ebnm_point_laplace)
    ) %>%
    flash_backfit(extrapolate = extrapolate, maxiter = 25, verbose = verbose)

  ### keep at most Kmax factors based on proportion of variance explained and refit EB-SNMF to gene expression data matrix
  kset <- (length(fit.snmf$pve) - rank(fit.snmf$pve) < Kmax) & (fit.snmf$pve > 0)
  fit.snmf <- flash_init(Y, S = 1/sqrt(nrow(Y)), var_type = 2) %>%
    flash_factors_init(
      init = lapply(fit.snmf$flash_fit$EF, function(x) x[, kset]),
      ebnm_fn = c(prior, ebnm_point_laplace)
    ) %>%
    flash_backfit(extrapolate = extrapolate, maxiter = maxiter, verbose = verbose)
}


### initialize the EB-SNMF fit from an EBMF fit without nonnegative constraints on L
init.ebsnmf <- function(fl, epsilon) {
  LL <- fl$flash_fit$EF[[1]]
  FF <- fl$flash_fit$EF[[2]]

  LL <- cbind(pmax(LL, 0), pmax(-LL, 0))
  FF <- cbind(FF, -FF)

  to.keep <- (colSums(LL) > epsilon)
  LL <- LL[, to.keep, drop = FALSE]
  FF <- FF[, to.keep, drop = FALSE]

  return(list(LL, FF))
}
