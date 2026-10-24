#!/usr/bin/env Rscript
library(MASS)
library(fda)
library(ClusterR)
library(flexmix)
library(Matrix)
library(aricode)

source("pre_clustering.R")
source("HFGM_ADMM.R")
source("HFLM_ADMM.R")
source("GHFM_testing.R")
source("Additional_Functions.R")

simulate_data <- function(response, n, K_true, beta_type, sigma_b, sigma_e, T_grid = 1440) {
  X_mat    <- matrix(rnorm(n * T_grid, 3, 1), n, T_grid)
  grp_true <- sample(1:K_true, n, replace = TRUE)
  bbasis   <- create.bspline.basis(c(0, T_grid), nbasis = 15, norder = 4)
  if (beta_type == "bspline") {
    mu <- c(5, -5, 2, -2)[1:K_true]
    if (response == "Bernoulli") {
      mu <- c(1.5, -1.5, 0.6, -0.6)[1:K_true]
    }
    beta_coefs <- matrix(rnorm(n * 15, 0, sigma_b), n, 15) + mu[grp_true]
  } else {
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
    beta_coefs <- t(coef_mat)[grp_true, , drop = FALSE]
  }
  fdX  <- Data2fd(y = t(X_mat), basisobj = bbasis,
                  argvals = seq(0, T_grid, length.out = T_grid))
  Zmat <- inprod(fdX$basis, bbasis)
  lin  <- rowSums((t(fdX$coefs) %*% Zmat) * beta_coefs) / T_grid
  if (response == "Gaussian") {
    Y <- lin + rnorm(n, 0, sigma_e)
  } else {
    Y <- rbinom(n, 1, 1 / (1 + exp(-lin)))
  }
  list(
    X_mat      = X_mat,
    fdX        = fdX,
    Y          = Y,
    grp_true   = grp_true,
    response   = response,
    true_basis = bbasis,
    beta_coef_mat_n_row = beta_coefs / T_grid
  )
}

fit_GHFM <- function(sim, K_pre, beta.basis, phi, eps0, lambda_rel, tau_q, theta, max_iter, C_n) {
  n <- length(sim$Y)
  L <- beta.basis$nbasis
  if (K_pre == 0) {
    pre_grp <- 1:n
  } else if (sim$response == "Gaussian") {
    pre_grp <- beta_clustering(K_pre, sim$X_mat, sim$Y, beta.basis, phi, eps0,
                               max.iteration = 50, tol = 0.001, x.basis.num = 15)
  } else {
    pre_grp <- beta_clustering_logistic(K_pre, sim$X_mat, sim$Y, beta.basis, phi, eps0,
                                        max.iteration = 50, tol = 0.001, x.basis.num = 15, max.iter = 25)
  }
  K <- length(unique(pre_grp))
  Gamma <- compute_Gamma_custom(sim$fdX, beta.basis)
  eta <- rep(0, n)
  w <- rep(1, n)
  z <- sim$Y
  for (it in 1:ifelse(sim$response == "Gaussian", 1, 25)) {
    if (sim$response == "Bernoulli") {
      p <- 1 / (1 + exp(-eta))
      w <- pmax(p * (1 - p), 0.001 * 0.999)
      z <- eta + (sim$Y - p) / w
    }
    b0 <- pre_clus_fit(Gamma, z, w, pre_grp, beta.basis, phi, eps0)$b.all
    eta <- rowSums(Gamma * b0[pre_grp, , drop = FALSE])
  }
  pairs <- combn(K, 2)
  dist0 <- sqrt(rowSums((b0[pairs[1, ], , drop = FALSE] - b0[pairs[2, ], , drop = FALSE])^2))
  best <- NULL
  for (tq in tau_q) {
    tau1 <- as.numeric(quantile(dist0, tq))
    for (lr in lambda_rel) {
      lam <- lr * max(dist0) * theta * tau1
      if (sim$response == "Gaussian") {
        admm <- hfgm_admm_with_pre_grp(sim$fdX, sim$Y, T = 1440, d = 3, M = L - 3, phi = phi / K, eps0 = eps0 / K,
                                       lambda = lam, theta = theta, tau1 = tau1, tau2 = 1e-4 * median(dist0),
                                       tau3 = 1e-2 * median(dist0), max_iter = max_iter, lambda0 = 0,
                                       init_b = as.vector(t(b0)), pre_grp = pre_grp)
      } else {
        admm <- hflm_admm_with_pre_clus(sim$fdX, sim$Y, T = 1440, d = 3, M = L - 3, phi = phi / K, eps0 = eps0 / K,
                                        lambda = lam, theta = theta, tau1 = tau1, tau2 = 1e-4 * median(dist0),
                                        tau3 = 1e-2 * median(dist0), max_iter = max_iter, lambda0 = 0,
                                        grp = pre_grp, init_b = as.vector(t(b0)))
      }
      final_grp <- admm$group[pre_grp]
      eta <- rowSums(Gamma * admm$b.group[final_grp, , drop = FALSE])
      if (sim$response == "Gaussian") {
        bic <- n * log(sum((sim$Y - eta)^2) / n) + C_n * L * length(unique(final_grp)) * log(n)
      } else {
        bic <- 2 * sum(log1p(exp(-abs(eta))) + pmax(eta, 0) - sim$Y * eta) + C_n * L * length(unique(final_grp)) * log(n)
      }
      if (is.null(best) || bic < best$bic) {
        best <- list(bic = bic, final_grp = final_grp, bhat_mat = admm$b.group[final_grp, , drop = FALSE])
      }
    }
  }
  list(pre_grp = pre_grp, final_grp = best$final_grp, bhat_mat = best$bhat_mat, Gamma = Gamma)
}

one_replication <- function(sc, rep, sc_id) {
  beta.basis <- create.bspline.basis(c(0, 1440), nbasis = 20, norder = 4)
  fm.basis   <- create.bspline.basis(c(0, 1440), nbasis = 5, norder = 4)
  sigma_b <- ifelse(sc$beta_type == "bspline", ifelse(sc$response == "Gaussian", 1, 0.3), 0)
  phi  <- 1e9
  eps0 <- ifelse(sc$response == "Gaussian", 1e5 / sc$n, 1e3)
  if (sc$n <= 500) {
    K_pre_list <- c(0, 10, 20)
    C_n <- 1
  } else {
    K_pre_list <- c(50, 100, 150)
    C_n <- log(sc$n)
  }
  set.seed(1000 * sc_id + rep)
  sim <- simulate_data(sc$response, sc$n, sc$K_true, sc$beta_type, sigma_b, 0.5)
  results <- list()
  for (K_pre in K_pre_list) {
    theta <- mean(colSums(compute_Gamma_custom(sim$fdX, beta.basis)^2)) / sc$n /
      ifelse(K_pre == 0, sc$n, K_pre) * ifelse(sc$response == "Gaussian", 1, 0.1)
    ghfm <- fit_GHFM(sim, K_pre, beta.basis, phi, eps0, lambda_rel = 10^seq(1.5, -1.5, by = -0.5),
                     tau_q = c(0.25, 0.5, 0.75), theta = theta, max_iter = 100, C_n = C_n)
    Gamma_test <- compute_Gamma_custom(sim$fdX, sim$fdX$basis)
    if (sc$response == "Gaussian") {
      P_val <- fit_full_vs_reduced_gauss(Gamma_test, sim$Y, ghfm$final_grp)$p_value
    } else {
      P_val <- suppressWarnings(fit_full_vs_reduced_ber(Gamma_test, sim$Y, ghfm$final_grp)$p_value)
    }
    results[[length(results) + 1]] <- data.frame(
      scenario_id = sc_id, response = sc$response, beta_type = sc$beta_type, n = sc$n, rep = rep,
      method  = ifelse(K_pre == 0, "GHFM*", paste0("GHFM_", K_pre)),
      NMI     = NMI(sim$grp_true, ghfm$final_grp, variant = "sum"),
      Purity  = ifelse(K_pre == 0, NA, compute_soft_purity(sim$grp_true, ghfm$pre_grp)),
      ISE     = compute_ISE(ghfm$bhat_mat, sim$beta_coef_mat_n_row, beta.basis, sim$true_basis),
      P_Value = P_val
    )
  }
  others <- list(
    Oracle = fit_Oracle(sim$fdX, sim$Y, beta.basis, phi, eps0, sc$response, sim$grp_true, 25, 1e-6),
    FunMR  = suppressWarnings(fit_FunMR(sim$fdX, sim$Y, fm.basis, sc$response, 6)),
    DevKM  = fit_DevKM(sim$fdX, sim$Y, beta.basis, phi, eps0, sc$response, 6, 10, 25, 1e-6),
    HSFLM  = fit_HSFLM(sim$fdX, sim$Y, beta.basis, phi, eps0, sc$response, 25, 1e-6)
  )
  for (m in names(others)) {
    results[[length(results) + 1]] <- data.frame(
      scenario_id = sc_id, response = sc$response, beta_type = sc$beta_type, n = sc$n, rep = rep,
      method  = m,
      NMI     = ifelse(m == "HSFLM", NA, NMI(sim$grp_true, others[[m]]$cls, variant = "sum")),
      Purity  = NA,
      ISE     = compute_ISE(others[[m]]$bhat_mat, sim$beta_coef_mat_n_row,
                            if (m == "FunMR") fm.basis else beta.basis, sim$true_basis),
      P_Value = NA
    )
  }
  do.call(rbind, results)
}

scenarios <- expand.grid(
  response  = c("Gaussian", "Bernoulli"),
  beta_type = c("bspline", "linear_quad"),
  n         = c(100, 10000),
  K_true    = 2,
  stringsAsFactors = FALSE
)
n_rep <- 1

results <- list()
sc_id <- 0
for (sc in split(scenarios, seq_len(nrow(scenarios)))) {
  sc_id <- sc_id + 1
  for (rep in 1:n_rep) {
    results[[length(results) + 1]] <- one_replication(sc, rep, sc_id)
  }
}

res_df <- do.call(rbind, results)
res_tab <- aggregate(cbind(NMI, Purity, P_Value) ~ scenario_id + response + beta_type + n + method,
                     data = res_df, FUN = function(x) round(mean(x), 3), na.action = na.pass)
res_tab$ISE <- aggregate(ISE ~ scenario_id + response + beta_type + n + method, data = res_df,
                         FUN = function(x) round(median(x), 3))$ISE
print(res_tab[order(res_tab$scenario_id, res_tab$method), ], row.names = FALSE)
