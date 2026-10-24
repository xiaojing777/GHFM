compute_ISE <- function(bhat_mat, btrue_mat, hat.basis, true.basis) {
  G_hh <- inprod(hat.basis, hat.basis)
  G_ht <- inprod(hat.basis, true.basis)
  G_tt <- inprod(true.basis, true.basis)
  err_sq   <- sum((bhat_mat %*% G_hh) * bhat_mat) - 2 * sum((bhat_mat %*% G_ht) * btrue_mat) +
    sum((btrue_mat %*% G_tt) * btrue_mat)
  true_sq  <- sum((btrue_mat %*% G_tt) * btrue_mat)
  sqrt(err_sq / true_sq)
}

compute_soft_purity <- function(true_labels, pred_labels) {
  clusters <- unique(pred_labels)
  p_k <- sapply(clusters, function(c) {
    idx <- pred_labels == c
    if (!any(idx)) return(0)
    tbl <- table(true_labels[idx])
    max(tbl) / sum(tbl)
  })
  mean(p_k)
}

fit_HSFLM <- function(X.fd, y, beta.basis, phi, eps0, response, max.iter, tol, intercept = FALSE) {
  if (response == "Gaussian") {
    fit <- ord_func_reg(X.fd, y, beta.basis, phi, eps0, intercept)
  } else {
    fit <- ord_func_logistic(X.fd, y, beta.basis, phi, eps0, max.iter, tol, intercept)
  }
  list(cls = rep(1, length(y)), alpha = fit$alpha,
       bhat_mat = matrix(as.vector(fit$b.hat), nrow = length(y), ncol = beta.basis$nbasis, byrow = TRUE))
}

fit_Oracle <- function(X.fd, y, beta.basis, phi, eps0, response, true_labels, max.iter, tol, intercept = FALSE) {
  bhat_mat <- matrix(0, nrow = length(y), ncol = beta.basis$nbasis)
  for (k in unique(true_labels)) {
    idk <- which(true_labels == k)
    bhat_mat[idk, ] <- fit_HSFLM(X.fd[idk], y[idk], beta.basis, phi, eps0, response, max.iter, tol, intercept)$bhat_mat
  }
  list(cls = true_labels, bhat_mat = bhat_mat)
}

fit_DevKM <- function(X.fd, y, beta.basis, phi, eps0, response, K_max, nstart, max.iter, tol, intercept = FALSE) {
  n <- length(y)
  homo <- fit_HSFLM(X.fd, y, beta.basis, phi, eps0, response, max.iter, tol, intercept)
  Gamma <- t(X.fd$coefs) %*% inprod(X.fd$basis, beta.basis)
  eta <- homo$alpha + rowSums(Gamma * homo$bhat_mat)
  if (response == "Gaussian") {
    dev <- y - eta
  } else {
    dev <- sign(y - 1 / (1 + exp(-eta))) * sqrt(2 * (log1p(exp(-abs(eta))) + pmax(eta, 0) - y * eta))
  }
  km_models <- lapply(1:K_max, function(k) kmeans(dev, centers = k, nstart = nstart))
  bic_vec <- sapply(km_models, function(km) {
    sig <- sqrt(km$tot.withinss / n)
    dens <- sapply(seq_along(km$size), function(k) km$size[k] / n * dnorm(dev, km$centers[k], sig))
    -2 * sum(log(rowSums(matrix(dens, nrow = n)))) + 2 * length(km$size) * log(n)
  })
  cls <- km_models[[which.min(bic_vec)]]$cluster
  bhat_mat <- matrix(0, nrow = n, ncol = beta.basis$nbasis)
  for (k in unique(cls)) {
    idk <- which(cls == k)
    bhat_mat[idk, ] <- fit_HSFLM(X.fd[idk], y[idk], beta.basis, phi, eps0, response, max.iter, tol, intercept)$bhat_mat
  }
  list(cls = cls, bhat_mat = bhat_mat)
}

fit_FunMR <- function(X.fd, y, beta.basis, response, K_max, intercept = FALSE) {
  Gamma <- t(X.fd$coefs) %*% inprod(X.fd$basis, beta.basis)
  L <- ncol(Gamma)
  df <- data.frame(y = y, Gamma)
  models <- lapply(1:K_max, function(k) tryCatch({
    if (response == "Gaussian") {
      flexmix(if (intercept) y ~ . else y ~ . - 1, data = df, k = k)
    } else {
      flexmix(if (intercept) cbind(y, 1 - y) ~ . else cbind(y, 1 - y) ~ . - 1, data = df, k = k,
              model = FLXMRglm(family = "binomial"))
    }
  }, error = function(e) NULL))
  models <- models[!sapply(models, is.null)]
  best <- models[[which.min(sapply(models, BIC))]]
  cls <- clusters(best)
  coefs <- parameters(best)
  bhat_mat <- t(coefs[paste0("coef.X", 1:L), cls, drop = FALSE])
  list(cls = cls, bhat_mat = bhat_mat)
}

make_beta_df <- function(coef_mat, bbasis, tag, scenario_id) {
  tgrid <- seq(0, 1440, length.out = 200)
  fdobj <- fd(coef_mat, bbasis)
  vals  <- eval.fd(tgrid, fdobj)
  K     <- ncol(vals)
  data.frame(
    time     = rep(tgrid, K),
    beta     = as.vector(vals),
    method   = tag,
    group_id = factor(rep(1:K, each = length(tgrid))),
    scenario = scenario_id
  )
}
