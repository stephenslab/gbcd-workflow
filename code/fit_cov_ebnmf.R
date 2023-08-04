### fit Empirical Bayes Matrix Factorization to single cell RNA-seq data based on covariance decomposition, with a nonnegative prior on GEP memberships
### Y: a cell by gene matrix of normalized and log-transformed gene expression data
### Kmax: a positive integer (at least 2) specifying an upper bound of the number of GEPs for initialization, which is approximately but often not exactly the final number of GEPs
### prior: a nonnegative prior for GEP memberships, such as generalized binary prior, which must be a function defined in the ebnm package
### thres: a positive number in (0,1) to specify threshold on Pearson correlation between L and L-tilde and to use concordant ones only
### extrapolate: a logical indicator specifying whether to use extrapolation to accelerate backfitting GEP memberships; see the flashier package for details
### maxiter: a positive integer specifying the maximum number of iterations to backfit GEP memberships
### verbose: an integer specifying whether and how to display progress updates, as described in the flashier package
fit.cov.ebnmf <- function(Y, Kmax, prior = ebnm::ebnm_generalized_binary, thres = 0.8, extrapolate = TRUE, maxiter = 500, verbose = 1){

  ### calculate the covariance matrix from the cell by gene matrix of gene expression data
  dat <- Y %*% t(Y)/ncol(Y)

  ### fit EBMF with point laplace prior to covariance matrix without considering the diagonal component
  fit.cov <- flash_init(dat, var_type = 0) %>%
    flash_greedy(Kmax = 1, ebnm_fn = ebnm_unimodal_nonnegative, verbose = 1) %>%
    flash_greedy(Kmax = Kmax - 1, ebnm_fn = ebnm_point_laplace, verbose = 1) %>%
    flash_backfit(maxiter = maxiter, verbose = 1)

  ### fit EBMF with point laplace prior to covariance matrix with the diagonal component
  fit.cov <- fit.ebcovmf(dat = dat, fl = fit.cov, prior = ebnm::ebnm_point_laplace, maxiter = maxiter, verbose = 1)$fl

  ### initialize EB-NMF fit from the EBMF fit with point laplace prior
  cov.init <- init.cov.ebnmf(fit.cov)
  kmax <- which.max(colSums(cov.init[[1]]))
  fit.cov <- flash_init(dat, var_type = 0) %>%
    flash_factors_init(
      init = lapply(cov.init, function(x) x[, kmax, drop = FALSE]),
      ebnm_fn = ebnm_unimodal_nonnegative
    ) %>%
    flash_factors_init(
      init = lapply(cov.init, function(x) x[, -c(kmax), drop = FALSE]),
      ebnm_fn = prior
    ) %>%
    flash_backfit(extrapolate = extrapolate, maxiter = 25, verbose = verbose)

  ### fit EB-NMF with a nonnegative prior to covariance matrix with the diagonal component
  fit.cov <- fit.ebcovmf(dat = dat, fl = fit.cov, prior = prior, extrapolate = extrapolate, maxiter = 25, verbose = verbose)$fl

  ### keep at most Kmax factors based on proportion of variance explained and refit EB-NMF to covariance matrix
  kset <- (length(fit.cov$pve) - rank(fit.cov$pve) < 1.5*Kmax) & (fit.cov$pve > 0)
  # kset <- (fit.cov$pve > 0)
  kall <- 1:fit.cov$n_factors
  if(!all(kset))
    fit.cov <- flash_factors_remove(fit.cov, kset=kall[!kset])
  fit.cov <- fit.ebcovmf(dat = dat, fl = fit.cov, prior = prior, extrapolate = extrapolate, maxiter = maxiter, verbose = verbose)$fl

  ### scale GEP membership estimates to 0-1 scale, and calculate Pearson correlations between L and L-tilde
  k.order <- order(fit.cov$pve, decreasing = TRUE)
  fit.L <- fit.cov$L_pm[, k.order]
  fit.Ltilde <- fit.cov$F_pm[, k.order]
  corr <- diag(cor(fit.L, fit.Ltilde))
  fit.L <- t(t(fit.L)/apply(fit.L, 2, max))

  ### estimate GEP signatures by fitting EB-SNMF to gene expression data matrix with fixed L estimated from covariance decomposition above
  init.F <- t(solve(crossprod(fit.L), crossprod(fit.L, Y)))
  fit.snmf <- flash_init(Y, S = 1/sqrt(nrow(Y)), var_type = 2) %>%
    flash_factors_init(
      init = list(as.matrix(fit.L), as.matrix(init.F)),
      ebnm_fn = c(prior, ebnm_point_laplace)
    ) %>%
    flash_factors_fix(kset = 1:ncol(fit.L), which_dim = "loadings") %>%
    flash_backfit(extrapolate = FALSE, verbose = verbose)

  ### calculate the z-score and lfsr of GEP signatures by running linear regression followed by ash
  J <- ncol(Y)
  genes <- colnames(Y)
  F.est <- matrix(0, J, ncol(fit.snmf$L_pm))
  rownames(F.est) <- genes
  colnames(F.est) <- colnames(fit.snmf$L_pm)
  F.se <- F.est
  F.z.ash <- F.est
  F.lfsr.ash <- F.est

  for (j in 1:J) {
    y     <- Y[,j]
    dat   <- as.data.frame(cbind(y, fit.snmf$L_pm))
    fit   <- lm(y ~ 0 + ., dat)
    coefs <- summary(fit)$coefficients
    F.est[j,] <- coefs[, "Estimate"]
    F.se[j,]  <- coefs[, "Std. Error"]
  }

  for(k in 1:ncol(F.est)){
    fit <- ash(F.est[,k], F.se[,k], mixcompdist = "normal", method = "shrink")
    F.z.ash[,k] <- fit$result$PosteriorMean/fit$result$PosteriorSD
    F.lfsr.ash[,k] <- fit$result$lfsr
  }

  ### return the estimated memberships and signatures only for GEPs whose L and Ltilde from covariance decomposition are highly concordant
  k.idx <- which(corr > thres)
  L.pm <- fit.snmf$L_pm[, k.idx]
  colnames(L.pm) <- c("Baseline", paste0("GEP", 1:(ncol(L.pm)-1)))
  F.lfc <- fit.snmf$F_pm[, k.idx]/log(2)
  F.z <- F.z.ash[, k.idx]
  F.lfsr <- F.lfsr.ash[, k.idx]
  colnames(F.lfc) <- colnames(L.pm)
  colnames(F.z) <- colnames(L.pm)
  colnames(F.lfsr) <- colnames(L.pm)
  return(list(L = L.pm, F = list(lfc = F.lfc, z_score = F.z, lfsr = F.lfsr)))
}


### initialize the EB-NMF fit to covariance matrix YY' s.t. E[YY'] = LL'+ D from an estimate of L without nonnegative constraints
init.cov.ebnmf <- function(fl, kset = 1:ncol(fl$flash_fit$EF[[1]])) {
  LL <- fl$flash_fit$EF[[1]][, kset, drop = FALSE]
  FF <- fl$flash_fit$EF[[2]][, kset, drop = FALSE]

  LL <- cbind(pmax(LL, 0), pmax(-LL, 0))
  LL <- cbind(fl$flash_fit$EF[[1]][, -kset, drop = FALSE], LL)
  FF <- cbind(pmax(FF, 0), pmax(-FF, 0))
  FF <- cbind(fl$flash_fit$EF[[2]][, -kset, drop = FALSE], FF)

  to.keep <- (colSums(LL) > .Machine$double.eps) & (colSums(FF) > .Machine$double.eps)
  LL <- LL[, to.keep, drop = FALSE]
  FF <- FF[, to.keep, drop = FALSE]

  return(list(LL, FF))
}


### fit EBMF to covariance matrix YY' s.t. E[YY'] = LL'+ D, where D = s2*I and I is an identity matrix
fit.ebcovmf <- function(dat, fl, prior, extrapolate = TRUE, maxiter = 500, epsilon = 0, verbose = 1){
  s2 <- max(0, mean(diag(dat) - diag(fitted(fl))))
  s2_diff <- Inf

  ### alternate between estimating s2 and backfitting until convergence
  while(s2 > 0 && abs(s2_diff - 1) > 1e-3) {
    dat_minuss2 <- dat - diag(rep(s2, ncol(dat)))
    kset <- which(fl$pve > epsilon)
    kmax <- which.max(fl$pve)
    kset <- kset[-c(kmax)]
    fl <- flash_init(dat_minuss2) %>%
      flash_factors_init(
        init = lapply(fl$flash_fit$EF, function(x) x[, kmax, drop = FALSE]),
        ebnm_fn = ebnm_unimodal_nonnegative
      ) %>%
      flash_factors_init(
        init = lapply(fl$flash_fit$EF, function(x) x[, kset, drop = FALSE]),
        ebnm_fn = prior
      ) %>%
      flash_backfit(extrapolate = extrapolate, maxiter = maxiter, verbose = verbose)
    old_s2 <- s2
    s2 <- max(0, mean(diag(dat) - diag(fitted(fl))))
    s2_diff <- s2 / old_s2
  }

  return(list(dat=dat, fl=fl, s2=s2))
}
