library(fda)
library(MASS)
library(ClusterR)

ord_func_reg = function(x.fd, y, beta.basis, phi, eps0, intercept = FALSE) {
  y <- matrix(y, ncol = 1)
  n <- nrow(y)
  Z0 <- inprod(x.fd$basis, beta.basis)
  Z <- t(x.fd$coefs) %*% Z0
  R <- eval.penalty(beta.basis, int2Lfd(2))
  Zc <- sweep(Z, 2, intercept * colMeans(Z))
  b.hat <- solve(t(Zc) %*% Zc / n + phi * R + eps0 * diag(ncol(Z))) %*% t(Zc) %*% (y - intercept * mean(y)) / n
  alpha <- intercept * mean(y - Z %*% b.hat)
  return(list(b.hat = b.hat, alpha = alpha))
}

pre_clus_fit = function(Gamma, z, w, grp, beta.basis, phi, eps0, intercept = FALSE) {
  R <- eval.penalty(beta.basis, int2Lfd(2))
  labs <- sort(unique(grp))
  b.d <- matrix(0, nrow = length(labs), ncol = ncol(Gamma))
  b.c <- b.d
  for (l in seq_along(labs)) {
    idx <- which(grp == labs[l])
    G.l <- Gamma[idx, , drop = FALSE]
    A.inv <- ginv(t(G.l) %*% (w[idx] * G.l) + length(idx) * (phi * R + eps0 * diag(ncol(Gamma))))
    b.d[l, ] <- A.inv %*% t(G.l) %*% (w[idx] * z[idx])
    b.c[l, ] <- A.inv %*% t(G.l) %*% w[idx]
  }
  k <- match(grp, labs)
  alpha <- 0
  if (intercept) {
    alpha <- (sum(w * z) - sum(w * rowSums(Gamma * b.d[k, , drop = FALSE]))) /
      (sum(w) - sum(w * rowSums(Gamma * b.c[k, , drop = FALSE])))
  }
  return(list(labs = labs, b.all = b.d - alpha * b.c, alpha = alpha))
}

beta_clustering = function(K, X, y, beta.basis, phi, eps0, max.iteration, tol, grp = NULL, x.basis.num,
                           intercept = FALSE) {
  n <- nrow(X)
  y <- as.numeric(y)
  x.basis <- create.bspline.basis(rangeval = c(0, ncol(X)), nbasis = x.basis.num, norder = 4)
  time_points <- seq(0, ncol(X), length.out = ncol(X))
  X.fd <- Data2fd(y = t(as.matrix(X)), basisobj = x.basis, argvals = time_points)
  X.fd.mat <- t(X.fd$coefs)
  if (is.null(grp)) {
    grp <- ClusterR::KMeans_rcpp(X.fd.mat, clusters = K, num_init = 5,
                                 max_iters = 100, initializer = "kmeans++")$clusters
  }
  Gamma <- X.fd.mat %*% inprod(X.fd$basis, beta.basis)
  it <- 1
  while (it <= max.iteration) {
    grp.old <- grp
    fit <- pre_clus_fit(Gamma, y, rep(1, n), grp, beta.basis, phi, eps0, intercept)
    distances <- (y - fit$alpha - Gamma %*% t(fit$b.all))^2
    grp <- fit$labs[apply(distances, 1, which.min)]
    if (sum(grp.old != grp) <= tol * n) break
    it <- it + 1
  }
  return(match(grp, sort(unique(grp))))
}

ord_func_logistic = function(x.fd, y, beta.basis, phi, eps0, max.iter, tol, intercept = FALSE) {
  y <- as.numeric(y)
  n <- length(y)
  Z0 <- inprod(x.fd$basis, beta.basis)
  gamma_mat <- t(x.fd$coefs) %*% Z0
  R <- eval.penalty(beta.basis, int2Lfd(2))
  b <- matrix(0, nrow = ncol(gamma_mat), ncol = 1)
  alpha <- 0
  eta <- rep(0, n)
  for (iter in 1:max.iter) {
    p <- 1 / (1 + exp(-eta))
    W <- pmax(p * (1 - p), 1e-4 * (1 - 1e-4))
    tilde_y <- eta + (y - p) / W
    gamma_bar <- intercept * colSums(W * gamma_mat) / sum(W)
    y_bar <- intercept * sum(W * tilde_y) / sum(W)
    gamma_c <- sweep(gamma_mat, 2, gamma_bar)
    lhs <- crossprod(gamma_c, W * gamma_c) / n + phi * R + eps0 * diag(ncol(gamma_mat))
    rhs <- crossprod(gamma_c, W * (tilde_y - y_bar)) / n
    b <- solve(lhs, rhs)
    alpha <- y_bar - sum(gamma_bar * b)
    eta_new <- as.vector(alpha + gamma_mat %*% b)
    if (max(abs(eta_new - eta)) < tol) break
    eta <- eta_new
  }
  return(list(b.hat = b, alpha = alpha))
}

beta_clustering_logistic = function(K_pre, X, y, beta.basis, phi, eps0, max.iteration, tol,
                                    grp_init = NULL, x.basis.num, max.iter, intercept = FALSE) {
  n <- nrow(X)
  y <- as.numeric(y)
  x.basis <- create.bspline.basis(rangeval = c(0, ncol(X)), nbasis = x.basis.num, norder = 4)
  time_points <- seq(0, ncol(X), length.out = ncol(X))
  X.fd <- Data2fd(y = t(as.matrix(X)), argvals = time_points, basisobj = x.basis)
  X.fd.mat <- t(X.fd$coefs)
  if (is.null(grp_init)) {
    grp_init <- ClusterR::KMeans_rcpp(X.fd.mat, clusters = K_pre, num_init = 5,
                                      max_iters = 100, initializer = "kmeans++")$clusters
  }
  Gamma <- X.fd.mat %*% inprod(X.fd$basis, beta.basis)
  current_grp <- grp_init
  eta <- rep(0, n)
  it <- 1
  while (it <= max.iteration) {
    for (iter in 1:max.iter) {
      p <- 1 / (1 + exp(-eta))
      W <- pmax(p * (1 - p), 1e-4 * (1 - 1e-4))
      fit <- pre_clus_fit(Gamma, eta + (y - p) / W, W, current_grp, beta.basis, phi, eps0, intercept)
      eta_new <- fit$alpha + rowSums(Gamma * fit$b.all[match(current_grp, fit$labs), , drop = FALSE])
      conv <- max(abs(eta_new - eta)) < tol
      eta <- eta_new
      if (conv) break
    }
    eta_all <- fit$alpha + Gamma %*% t(fit$b.all)
    log_lik <- log1p(exp(-abs(eta_all))) + pmax(eta_all, 0) - y * eta_all
    new_grp <- fit$labs[apply(log_lik, 1, which.min)]
    change_num <- sum(new_grp != current_grp)
    current_grp <- new_grp
    eta <- fit$alpha + rowSums(Gamma * fit$b.all[match(current_grp, fit$labs), , drop = FALSE])
    if (change_num <= tol * n) break
    it <- it + 1
  }
  return(match(current_grp, sort(unique(current_grp))))
}
