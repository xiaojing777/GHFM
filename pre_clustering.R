library(fda)
library(MASS)
library(ClusterR)

ord_func_reg = function(x.fd, y, beta.basis, phi) {
  y <- matrix(y, ncol = 1)
  Z0 <- inprod(x.fd$basis, beta.basis)
  Z <- t(x.fd$coefs) %*% Z0
  R <- eval.penalty(beta.basis, int2Lfd(2))
  diag_adjust <- 0.01 * diag(nrow(R))
  b.hat <- solve(t(Z) %*% Z + phi * R + diag_adjust) %*% t(Z) %*% y
  X.train.26 <- t(x.fd$coefs)
  alpha <- mean(y - X.train.26 %*% Z0 %*% b.hat)
  return(list(b.hat = b.hat, alpha = alpha))
}

beta_clustering = function(K = 100, X, y, beta.basis, phi, max.iteration = 200,
                           tol = 0.005, grp = NULL, x.basis.num = 15) {
  n <- nrow(X)
  y <- matrix(y, ncol = 1)
  if (is.null(grp)) {
    k2 <- ClusterR::KMeans_rcpp(X, clusters = K, num_init = 5,
                                max_iters = 100, initializer = "kmeans++")
    grp <- k2$clusters
  }
  x.basis <- create.bspline.basis(rangeval = c(0, ncol(X)),
                                  nbasis = x.basis.num, norder = 4)
  time_points <- seq(0, ncol(X), length.out = ncol(X))
  X.mat <- as.matrix(X)
  X.fd <- Data2fd(y = t(X.mat), basisobj = x.basis, argvals = time_points)
  X.fd.mat <- t(X.fd$coefs)
  Z0 <- inprod(X.fd$basis, beta.basis)
  it <- 1
  while (it < max.iteration) {
    grp.old <- grp
    b.all <- matrix(NA, nrow = K, ncol = beta.basis$nbasis)
    for (l in 1:K) {
      idx <- which(grp == l)
      if (length(idx) == 0) next
      X.mat.l <- t(X.mat[idx, , drop = FALSE])
      X.fd.l <- Data2fd(y = X.mat.l, basisobj = x.basis, argvals = time_points)
      y.l <- matrix(y[idx], ncol = 1)
      reg_result <- ord_func_reg(X.fd.l, y.l, beta.basis, phi)
      b.all[l, ] <- reg_result$b.hat
    }
    valid_clusters <- which(!is.na(b.all[, 1]))
    b.all <- b.all[valid_clusters, , drop = FALSE]
    for (i in 1:n) {
      distances <- sapply(1:nrow(b.all), function(j) {
        x_vec <- matrix(X.fd.mat[i, ], nrow = 1)
        b_vec <- matrix(b.all[j, ], ncol = 1)
        tmp <- x_vec %*% Z0 %*% b_vec
        (y[i, 1] - tmp[1, 1])^2
      })
      grp[i] <- valid_clusters[which.min(distances)]
    }
    if (it >= 35 && sum(grp.old != grp) <= tol * n) break
    it <- it + 1
  }
  return(grp)
}

ord_func_logistic = function(x.fd, y, beta.basis, phi, max.iter = 50, tol = 1e-3) {
  y <- as.numeric(y)
  n <- length(y)
  Z0 <- inprod(x.fd$basis, beta.basis)
  gamma_mat <- t(x.fd$coefs) %*% Z0
  R <- eval.penalty(beta.basis, int2Lfd(2))
  diag_adjust <- 1e-3 * diag(nrow(R))
  b <- matrix(0, nrow = ncol(gamma_mat), ncol = 1)
  alpha <- 0
  for (iter in 1:max.iter) {
    eta <- as.vector(alpha + gamma_mat %*% b)
    p <- 1 / (1 + exp(-eta))
    p <- pmin(pmax(p, 1e-4), 1 - 1e-4)
    W <- p * (1 - p)
    W_sqrt <- sqrt(W)
    tilde_y <- eta + (y - p) / W
    bar_Y <- W_sqrt * tilde_y
    bar_gamma <- gamma_mat * W_sqrt
    lhs <- crossprod(bar_gamma) / n + phi * R + diag_adjust
    rhs <- crossprod(bar_gamma, bar_Y) / n
    b_new <- solve(lhs, rhs)
    alpha_new <- mean(bar_Y - bar_gamma %*% b_new)
    if (max(abs(b_new - b)) < tol && abs(alpha_new - alpha) < tol) {
      b <- b_new
      alpha <- alpha_new
      break
    }
    b <- b_new
    alpha <- alpha_new
  }
  return(list(b.hat = b, alpha = alpha))
}

beta_clustering_logistic = function(K_pre = 30, X, y, beta.basis, phi, max.iteration = 100,
                                    tol = 0.005, grp_init = NULL, x.basis.num = 15) {
  n <- nrow(X)
  y <- as.numeric(y)
  if (is.null(grp_init)) {
    grp_init <- ClusterR::KMeans_rcpp(X, clusters = K_pre)$clusters
  }
  x.basis <- create.bspline.basis(rangeval = c(0, ncol(X)), nbasis = x.basis.num, norder = 4)
  time_points <- seq(0, ncol(X), length.out = ncol(X))
  X.fd <- Data2fd(y = t(X), argvals = time_points, basisobj = x.basis)
  Z0 <- inprod(X.fd$basis, beta.basis)
  X.fd.mat <- t(X.fd$coefs)
  it <- 1
  current_grp <- grp_init
  new_grp <- current_grp
  while (it <= max.iteration) {
    b.all <- matrix(NA, nrow = K_pre, ncol = beta.basis$nbasis)
    alpha.all <- rep(NA, K_pre)
    for (k in 1:K_pre) {
      idx <- which(current_grp == k)
      if (length(idx) < 1) {
        b.all[k, ] <- NA_real_
        alpha.all[k] <- NA_real_
        next
      }
      X.l <- t(X[idx, , drop = FALSE])
      X.fd.l <- Data2fd(y = X.l, argvals = time_points, basisobj = x.basis)
      fit <- ord_func_logistic(X.fd.l, y[idx], beta.basis, phi = phi)
      b.all[k, ] <- fit$b.hat
      alpha.all[k] <- fit$alpha
    }
    valid_clusters <- which(!is.na(b.all[, 1]))
    if (length(valid_clusters) == 0) stop("all clusters failed")
    b.valid <- b.all[valid_clusters, , drop = FALSE]
    alpha.valid <- alpha.all[valid_clusters]
    for (i in 1:n) {
      gamma_i <- X.fd.mat[i, ] %*% Z0
      eta_all <- gamma_i %*% t(b.valid) + matrix(alpha.valid, nrow = 1)
      p_all <- 1 / (1 + exp(-eta_all))
      p_all <- pmin(pmax(p_all, 1e-3), 1 - 1e-3)
      log_lik <- - (y[i] * log(p_all) + (1 - y[i]) * log(1 - p_all))
      idx_increa_order <- order(log_lik)
      worst_idx        <- tail(idx_increa_order, 2)
      if (!(current_grp[i] %in% worst_idx)) {
        new_grp[i] <- current_grp[i]
        next
      }
      new_grp[i] <- valid_clusters[which.min(log_lik)]
    }
    change_num <- sum(new_grp != current_grp)
    current_grp <- new_grp
    it <- it + 1
  }
  return(current_grp)
}
