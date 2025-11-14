#!/usr/bin/env Rscript
library(MASS)
library(fda)
library(ClusterR)
library(flexmix)
library(mclust)
library(Matrix)
library(ggplot2)
library(dplyr)
library(tidyr)
library(aricode)

source("/Users/xiaojing/Desktop/GitHub Projects/GHFM/pre_clustering.R")
source("/Users/xiaojing/Desktop/GitHub Projects/GHFM/HFGM_ADMM.R")
source("/Users/xiaojing/Desktop/GitHub Projects/GHFM/HFLM_ADMM.R")
source("/Users/xiaojing/Desktop/GitHub Projects/GHFM/GHFM_testing.R")
source("/Users/xiaojing/Desktop/GitHub Projects/GHFM/Additional_Functions.R")



simulate_data <- function(response, n, K_true, K_pre, beta_type, sigma, T_grid = 1440) {
  X_mat    <- matrix(rnorm(n*T_grid,3,1), n, T_grid)
  grp_true <- rep(1:K_true, each = n/K_true)
  if (K_pre == 0) {
    init_grp <- 1:n
  } else {
    init_grp <- rep(1:K_pre, each = n/K_pre)
  }
  if (beta_type == "bspline") {
    bbasis   <- create.bspline.basis(c(0,T_grid), nbasis=15, norder=4)
    set.seed(100+sigma)
    coef_mat <- matrix(0,15,K_true)
    coef_mat[,1] <- rnorm(15,  5, sigma)
    coef_mat[,2] <- rnorm(15, -5, sigma)
  } else {
    bbasis   <- create.bspline.basis(c(0,T_grid), nbasis=15, norder=4)
    arg_long <- seq(0, T_grid, length.out = T_grid)
    coef_mat <- matrix(NA, nrow = 15, ncol = K_true)
    for (k in seq_len(K_true)) {
      if (k %% 2 == 1) {
        idx_lin     <- (k + 1) / 2
        intercept_k <- -3 + 0.5 * (idx_lin - 1)
        slope_k     <- 0.005 - 0.002 * (idx_lin - 1)
        f_k         <- intercept_k + slope_k * arg_long
      } else {
        idx_quad <- k / 2
        A <- -1 - 0.5  * (idx_quad - 1)
        B <- -0.005 + 0.001 * (idx_quad - 1)
        C <- 0.000005 - 0.000002 * (idx_quad - 1)
        f_k       <- A + B * arg_long + C * arg_long^2
      }
      fd_k         <- Data2fd(y = f_k, argvals = arg_long, basisobj = bbasis)
      coef_mat[, k] <- as.numeric(fd_k$coefs)
    }
  }
  p <- nrow(coef_mat)
  beta_coefs <- matrix(0, nrow = n, ncol = p)
  for (k in seq_len(K_true)) {
    rows <- which(grp_true == k)
    beta_coefs[rows, ] <- matrix(
      coef_mat[, k],
      nrow  = length(rows),
      ncol  = p,
      byrow = TRUE
    )
  }
  fdX  <- Data2fd(y=t(X_mat), basisobj=bbasis,
                  argvals=seq(0,T_grid,length.out=T_grid))
  Zmat <- inprod(fdX$basis, bbasis)
  if (response=="Gaussian") {
    lin <- rowSums((t(fdX$coefs)%*%Zmat)*beta_coefs)
    Y   <- lin + rnorm(n)
  } else {
    linp <- rowSums((t(fdX$coefs)%*%Zmat)*beta_coefs)-10
    p    <- 1/(1+exp(-linp))
    Y    <- rbinom(n,1,p)
  }
  list(
    X_mat    = X_mat,
    fdX      = fdX,
    Zmat     = Zmat,
    Y        = Y,
    grp_true = grp_true,
    init_grp = init_grp,
    coef_mat = coef_mat,
    response = response,
    beta_coef_mat_n_row = beta_coefs
  )
}

scenarios <- expand.grid(
  response  = c("Gaussian"),
  n         = c(100),
  K_true    = 2,
  K_pre     = c(10, 20),
  sigma     = c(0.1),
  beta_type = c("bspline"),
  stringsAsFactors = FALSE
)

results <- list()
beta_curves <- list()
sc_id <- 0
if (!exists("coef_storage")) coef_storage <- list()

for (sc in split(scenarios, seq_len(nrow(scenarios)))) {
  sc_id <- sc_id + 1
  sim <- simulate_data(sc$response, sc$n, sc$K_true, sc$K_pre, sc$beta_type, sc$sigma)
  if (sc$K_pre == 0) {
    pre_grp <- 1:sc$n
  } else if (sc$response == "Gaussian") {
    pre_grp <- beta_clustering(sc$K_pre, sim$X_mat, sim$Y, sim$fdX$basis,
                               phi=1, max.iteration=50, tol=0.005,
                               grp=sim$init_grp)
  } else {
    pre_grp <- beta_clustering_logistic(K_pre=sc$K_pre,
                                        X=sim$X_mat, y=sim$Y,
                                        beta.basis=sim$fdX$basis,
                                        phi=1, x.basis.num=15,
                                        grp_init=sim$init_grp)
  }
  dd <- t(sim$coef_mat)
  half <- length(unique(pre_grp))/2
  dd_new <- rbind(
    matrix(rep(dd[1, ], each = half), nrow = half, byrow = TRUE),
    matrix(rep(dd[2, ], each = half), nrow = half, byrow = TRUE)
  )
  dd_new_vec <- matrix(as.vector(t(dd_new)), ncol = 1)
  if (sc$response == "Bernoulli") {
    admm <- hflm_admm_with_pre_clus(
      X.fd      = sim$fdX, Y      = sim$Y,
      T         = 1440,           d = 3, M = 12,
      phi       = 1,             lambda = 0.001,
      theta     = 1,             tau1   = 1,
      tau2      = 1,             tau3   = 1,
      max_iter  = 50,            grp    = pre_grp
    )
  }
  if (sc$response == "Gaussian") {
    admm <- hfgm_admm_with_pre_grp(
      X.fd      = sim$fdX, Y      = sim$Y,
      T         = 1440,           d = 3, M = 12,
      phi       = 1,             lambda = 0.001,
      theta     = 1,             tau1   = 1,
      tau2      = 1,
      max_iter  = 50,            pre_grp    = pre_grp
    )
  }
  b_vec <- admm$b
  K     <- length(unique(pre_grp))
  L     <- sim$fdX$basis$nbasis
  bmat  <- matrix(b_vec,
                  nrow   = K,
                  ncol   = L,
                  byrow  = TRUE)
  if (sc$response=="Bernoulli") {
    fuse  <- fuse_group(sim$fdX$basis, bmat, tol = 0.03)
  } else {
    fuse  <- fuse_group(sim$fdX$basis, bmat, tol = 40)
  }
  final_grp_hflm <- fuse$group[pre_grp]
  K_est_hflm     <- length(unique(final_grp_hflm))
  coef_storage[[paste0("scenario_", sc_id)]] <- list(
    b.hat         = fuse$b.hat,
    true_coef_mat = sim$coef_mat,
    bmat          = bmat
  )
  G <- as.matrix(inprod(sim$fdX$basis, sim$fdX$basis))
  btrue_mat <- sim$beta_coef_mat_n_row
  cls <- final_grp_hflm
  bhat_mat <- fuse$b.hat[cls, , drop = FALSE]
  nmii   <- NMI(sim$grp_true, cls)
  sp     <- compute_soft_purity(sim$grp_true, cls)
  err    <- compute_ISE(bhat_mat, btrue_mat, G)
  if (sc$response == "Gaussian") {
    P_val <- fit_full_vs_reduced_gauss(t(sim$fdX$coefs)%*%sim$Zmat, sim$Y, cls)$p_value
  } else {
    P_val <- fit_full_vs_reduced_ber(t(sim$fdX$coefs)%*%sim$Zmat, sim$Y, cls)$p_value
  }
  results[[length(results) + 1]] <- data.frame(
    scenario_id = sc_id,
    response    = sc$response,
    n           = sc$n,
    K_true      = sc$K_true,
    K_pre       = sc$K_pre,
    sigma       = sc$sigma,
    beta_type   = sc$beta_type,
    method      = "GHFM",
    NMI         = nmii,
    Purity  = sp,
    ISE_ratio   = err,
    P_Value     = P_val
  )
}

res_df <- do.call(rbind, results)
print(res_df)
