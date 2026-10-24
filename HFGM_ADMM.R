library(fda)
library(Matrix)

compute_basis_custom <- function(T, d, M) {
  create.bspline.basis(rangeval = c(0, T), nbasis = M + d, norder = d + 1)
}

compute_R_custom <- function(basis, w) {
  R1 <- eval.penalty(basis, int2Lfd(2))
  array(sapply(w, function(m) m * R1), dim = c(nrow(R1), ncol(R1), length(w)))
}

compute_Gamma_custom <- function(X.fd, basis) {
  X.coef <- t(X.fd$coefs)
  inprod_mat <- inprod(X.fd$basis, basis)
  Gamma <- X.coef %*% inprod_mat
  Gamma
}

compute_H_block_with_pre_grps <- function(Gamma, pre_grp) {
  n <- nrow(Gamma)
  L <- ncol(Gamma)
  k <- match(pre_grp, sort(unique(pre_grp)))
  sparseMatrix(i = rep(1:n, each = L), j = rep((k - 1) * L, each = L) + rep(1:L, n),
               x = as.vector(t(Gamma)), dims = c(n, max(k) * L))
}

compute_Omega_custom <- function(n, L) {
  cmat <- combn(n, 2)
  nrow_epsilon <- ncol(cmat)
  row_inds <- rep(seq_len(nrow_epsilon), each = 2)
  col_inds <- as.vector(rbind(cmat[1, ], cmat[2, ]))
  vals <- as.vector(rbind(rep(1, nrow_epsilon), rep(-1, nrow_epsilon)))
  eps_sparse <- sparseMatrix(i = row_inds, j = col_inds, x = vals, dims = c(nrow_epsilon, n))
  kronecker(eps_sparse, Diagonal(L))
}

shrinkage <- function(w, t) {
  nw <- sqrt(colSums(w^2))
  w * rep(ifelse(nw > t, 1 - t / nw, 0), each = nrow(w))
}

admm_factor <- function(B, theta) {
  L <- dim(B)[1]
  K <- dim(B)[3]
  Dinv <- vector("list", K)
  DB <- matrix(0, L, L)
  for (k in 1:K) {
    Dinv[[k]] <- solve(B[, , k] + theta * K * diag(L))
    DB <- DB + Dinv[[k]] %*% B[, , k]
  }
  list(Dinv = bdiag(Dinv), C = theta * K * solve(DB), L = L, K = K)
}

admm_solve <- function(fac, rhs, a, s) {
  Minv <- function(y) {
    z <- as.vector(fac$Dinv %*% y)
    u <- as.vector(fac$C %*% rowSums(matrix(z, nrow = fac$L)))
    z + as.vector(fac$Dinv %*% rep(u, fac$K))
  }
  if (length(rhs) == fac$K * fac$L) {
    return(Minv(rhs))
  }
  My <- Minv(rhs[-1])
  Ms <- Minv(s)
  x0 <- (rhs[1] - sum(s * My)) / (a - sum(s * Ms))
  c(x0, My - x0 * Ms)
}

hfgm_admm_without_pre_grp <- function(X.fd, Y, T, d, M, phi, eps0, lambda, theta, tau1, tau2, tau3, max_iter,
                                      lambda0, init_b = NULL, intercept = FALSE) {
  hfgm_admm_with_pre_grp(X.fd, Y, T, d, M, phi, eps0, lambda, theta, tau1, tau2, tau3, max_iter,
                         lambda0, init_b = init_b, pre_grp = 1:length(Y), intercept = intercept)
}

hfgm_admm_with_pre_grp <- function(X.fd, Y, T, d, M, phi, eps0, lambda, theta, tau1, tau2, tau3, max_iter,
                                   lambda0, init_b = NULL, pre_grp, intercept = FALSE) {
  Y <- matrix(Y, ncol = 1)
  n <- nrow(Y)
  basis <- compute_basis_custom(T, d, M)
  L <- basis$nbasis
  labs <- sort(unique(pre_grp))
  K <- length(labs)
  kidx <- match(pre_grp, labs)
  nk <- as.numeric(table(factor(pre_grp, levels = labs)))
  R <- compute_R_custom(basis, rep(1, K))
  Gamma <- compute_Gamma_custom(X.fd, basis)
  H <- compute_H_block_with_pre_grps(Gamma, pre_grp)
  Omega <- compute_Omega_custom(K, L)
  pairs <- combn(K, 2)
  B <- R
  for (k in 1:K) {
    B[, , k] <- crossprod(Gamma[kidx == k, , drop = FALSE]) / n + phi * R[, , k] + eps0 * diag(L)
  }
  s_vec <- as.vector(sapply(1:K, function(k) colSums(Gamma[kidx == k, , drop = FALSE]))) / n
  if (intercept) {
    H <- cbind(1, H)
    Omega <- cbind(0, Omega)
  }
  HtY <- as.vector((1 / n) * crossprod(H, Y))
  if (is.null(init_b)) {
    b <- admm_solve(admm_factor(B, lambda0), HtY, 1, s_vec)
  } else {
    b <- as.vector(init_b)
  }
  fac <- admm_factor(B, theta)
  alpha <- matrix(0, nrow = L, ncol = ncol(pairs))
  v <- matrix(0, nrow = L, ncol = ncol(pairs))
  for (m in 1:max_iter) {
    alpha_old <- alpha
    diff_b <- matrix(as.vector(Omega %*% b), nrow = L)
    thresh <- lambda / (theta * tau1) * (sqrt(colSums(alpha_old^2)) <= tau1)
    alpha <- shrinkage(diff_b + v / theta, thresh)
    b <- admm_solve(fac, HtY + as.vector(crossprod(Omega, as.vector(theta * alpha - v))), 1, s_vec)
    diff_b <- matrix(as.vector(Omega %*% b), nrow = L)
    v <- v + theta * (diff_b - alpha)
    r <- sqrt(sum((diff_b - alpha)^2))
    s <- theta * sqrt(sum(crossprod(Omega, as.vector(alpha - alpha_old))^2))
    if (r < tau2 && s < tau2) {
      break
    }
  }
  b0 <- 0
  if (intercept) {
    b0 <- b[1]
    b <- b[-1]
  }
  D <- matrix(0, K, K)
  D[t(pairs)] <- sqrt(colSums(alpha^2))
  group <- cutree(hclust(as.dist(D + t(D)), method = "single"), h = tau3)
  b_mat <- matrix(b, nrow = K, byrow = TRUE)
  b_group <- t(sapply(sort(unique(group)), function(g)
    colSums(nk[group == g] * b_mat[group == g, , drop = FALSE]) / sum(nk[group == g])))
  list(b = b, intercept = b0, alpha = as.vector(alpha), v = as.vector(v), basis = basis,
       group = group, b.group = b_group, iter = m)
}

fuse_group <- function(beta.basis, b.mat, tol) {
  Z <- inprod(beta.basis, beta.basis)
  C <- chol(Z)
  grp <- cutree(hclust(dist(b.mat %*% t(C)), method = "single"), h = tol)
  for (g in unique(grp)) {
    b.mat[grp == g, ] <- matrix(colMeans(b.mat[grp == g, , drop = FALSE]), nrow = sum(grp == g),
                                ncol = ncol(b.mat), byrow = TRUE)
  }
  list(group = grp, b.hat = b.mat)
}
