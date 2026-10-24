library(fda)
library(Matrix)

compute_basis_custom <- function(T, d, M) {
  create.bspline.basis(rangeval = c(0, T), nbasis = M + d, norder = d + 1)
}

compute_R_custom <- function(basis, w) {
  R1 <- eval.penalty(basis, int2Lfd(2))
  array(sapply(w, function(m) m * R1), dim = c(nrow(R1), ncol(R1), length(w)))
}

compute_Gamma_with_pre_clus_custom <- function(X.fd, basis, grp) {
  n <- ncol(X.fd$coefs)
  L <- basis$nbasis
  X.coef <- t(X.fd$coefs)
  inprod_mat <- inprod(X.fd$basis, basis)
  X.neiji <- X.coef %*% inprod_mat
  k <- match(grp, sort(unique(grp)))
  sparseMatrix(i = rep(1:n, each = L), j = rep((k - 1) * L, each = L) + rep(1:L, n),
               x = as.vector(t(X.neiji)), dims = c(n, max(k) * L))
}

compute_Omega_custom <- function(K, L) {
  cmat <- combn(K, 2)
  m <- ncol(cmat)
  row_inds <- rep(seq_len(m), each = 2)
  col_inds <- as.vector(rbind(cmat[1, ], cmat[2, ]))
  vals <- as.vector(rbind(rep(1, m), rep(-1, m)))
  eps_sparse <- sparseMatrix(i = row_inds, j = col_inds, x = vals, dims = c(m, K))
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

hflm_admm_with_pre_clus <- function(X.fd, Y, T, d, M, phi, eps0, lambda, theta,
                                    tau1, tau2, tau3, max_iter, lambda0, grp, init_b = NULL, intercept = FALSE) {
  Y <- as.vector(Y)
  n <- length(Y)
  labs <- sort(unique(grp))
  K <- length(labs)
  kidx <- match(grp, labs)
  nk <- as.numeric(table(factor(grp, levels = labs)))
  basis <- compute_basis_custom(T, d, M)
  L <- basis$nbasis
  R <- compute_R_custom(basis, rep(1, K))
  Gamma <- compute_Gamma_with_pre_clus_custom(X.fd, basis, grp)
  X.neiji <- t(X.fd$coefs) %*% inprod(X.fd$basis, basis)
  Omega <- compute_Omega_custom(K, L)
  pairs <- combn(K, 2)
  if (intercept) {
    Gamma <- cbind(1, Gamma)
    Omega <- cbind(0, Omega)
  }
  irls_blocks <- function(w, z) {
    B <- R
    for (k in 1:K) {
      idx <- which(kidx == k)
      B[, , k] <- crossprod(X.neiji[idx, , drop = FALSE] * sqrt(w[idx])) / n + phi * R[, , k] + eps0 * diag(L)
    }
    s_vec <- as.vector(sapply(1:K, function(k) colSums(X.neiji[kidx == k, , drop = FALSE] * w[kidx == k]))) / n
    list(B = B, s = s_vec, a = sum(w) / n, rhs = as.vector(crossprod(Gamma, w * z)) / n)
  }
  if (is.null(init_b)) {
    blk <- irls_blocks(rep(1, n), Y)
    b <- admm_solve(admm_factor(blk$B, lambda0), blk$rhs, blk$a, blk$s)
  } else {
    b <- as.vector(init_b)
  }
  alpha <- matrix(0, nrow = L, ncol = ncol(pairs))
  v <- matrix(0, nrow = L, ncol = ncol(pairs))
  for (m in 1:max_iter) {
    eta <- as.vector(Gamma %*% b)
    p <- 1 / (1 + exp(-eta))
    w <- pmax(p * (1 - p), 0.001 * 0.999)
    tilde_y <- eta + (Y - p) / w
    blk <- irls_blocks(w, tilde_y)
    alpha_old <- alpha
    diff_b <- matrix(as.vector(Omega %*% b), nrow = L)
    thresh <- lambda / (theta * tau1) * (sqrt(colSums(alpha_old^2)) <= tau1)
    alpha <- shrinkage(diff_b + v / theta, thresh)
    b <- admm_solve(admm_factor(blk$B, theta), blk$rhs + as.vector(crossprod(Omega, as.vector(theta * alpha - v))),
                    blk$a, blk$s)
    diff_b <- matrix(as.vector(Omega %*% b), nrow = L)
    v <- v + theta * (diff_b - alpha)
    r <- sqrt(sum((diff_b - alpha)^2))
    s <- theta * sqrt(sum(crossprod(Omega, as.vector(alpha - alpha_old))^2))
    if (r < tau2 && s < tau2) break
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
