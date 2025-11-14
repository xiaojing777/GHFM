library(fda)
library(Matrix)

compute_basis_custom <- function(T, d, M) {
  create.bspline.basis(rangeval = c(0, T), nbasis = M + d, norder = d + 1)
}

compute_R_custom <- function(basis, K) {
  R1 <- eval.penalty(basis, int2Lfd(2))
  rlist <- replicate(K, R1, simplify = FALSE)
  R_block <- bdiag(rlist)
  as.matrix(R_block)
}

compute_Gamma_with_pre_clus_custom <- function(X.fd, basis, grp) {
  n <- ncol(X.fd$coefs)
  L <- basis$nbasis
  X.coef <- t(X.fd$coefs)
  inprod_mat <- inprod(X.fd$basis, basis)
  X.neiji <- X.coef %*% inprod_mat
  K <- length(unique(grp))
  Gamma <- matrix(0, nrow = n, ncol = K * L)
  for (k in 1:K) {
    idx_k <- which(grp == k)
    if (length(idx_k) == 0) next
    neiji.k.data <- X.neiji[idx_k, , drop = FALSE]
    for (i in seq_along(idx_k)) {
      Gamma[idx_k[i], ((k - 1) * L + 1):(k * L)] <- neiji.k.data[i, ]
    }
  }
  Gamma
}

compute_Omega_custom <- function(K, L) {
  cmat <- combn(K, 2)
  m <- ncol(cmat)
  row_inds <- rep(seq_len(m), each = 2)
  col_inds <- as.vector(rbind(cmat[1, ], cmat[2, ]))
  vals <- as.vector(rbind(rep(1, m), rep(-1, m)))
  eps_sparse <- sparseMatrix(i = row_inds, j = col_inds, x = vals, dims = c(m, K))
  kronecker(eps_sparse, diag(L))
}

shrinkage <- function(w, t) {
  nw <- sqrt(sum(w^2))
  if (nw > t) {
    (1 - t / nw) * w
  } else {
    rep(0, length(w))
  }
}

hflm_admm_with_pre_clus <- function(X.fd, Y, T, d, M, phi, lambda, theta,
                                    tau1 = 1, tau2 = 1, tau3 = 1,
                                    max_iter = 100, grp, init_b = NULL) {
  Y <- matrix(Y, ncol = 1)
  n <- nrow(Y)
  K <- length(unique(grp))
  basis <- compute_basis_custom(T, d, M)
  L <- basis$nbasis
  R <- compute_R_custom(basis, K)
  Gamma <- compute_Gamma_with_pre_clus_custom(X.fd, basis, grp)
  Omega <- compute_Omega_custom(K, L)
  matrix_to_invert <- (1 / K) * crossprod(Gamma) + phi * R + lambda * crossprod(Omega)
  diag(matrix_to_invert) <- diag(matrix_to_invert) + 1e-3
  inv_matrix <- chol2inv(chol(matrix_to_invert))
  if (is.null(init_b)) {
    b <- inv_matrix %*% ((1 / K) * t(Gamma) %*% Y)
  } else {
    b <- init_b
  }
  b <- as.matrix(b)
  alpha <- matrix(0, nrow = K * (K - 1) / 2, ncol = L)
  v <- matrix(0, nrow = K * (K - 1) / 2, ncol = L)
  for (m in 1:max_iter) {
    p <- 1 / (1 + exp(-Gamma %*% b))
    p[p > 0.999] <- 0.999
    p[p < 0.001] <- 0.001
    w <- as.vector(p * (1 - p))
    bar_Gamma <- Gamma * sqrt(w)
    tilde_y <- Gamma %*% b + (Y - p) / w
    bar_Y <- sqrt(w) * tilde_y
    idxcount <- 0
    for (i in 1:(K - 1)) {
      for (j in (i + 1):K) {
        idxcount <- idxcount + 1
        b_i <- b[((i - 1) * L + 1):(i * L), , drop = FALSE]
        b_j <- b[((j - 1) * L + 1):(j * L), , drop = FALSE]
        diff_b <- as.vector(b_i - b_j)
        alpha[idxcount, ] <- shrinkage(diff_b + v[idxcount, ] / theta, lambda / theta)
      }
    }
    idxcount <- 0
    for (i in 1:(K - 1)) {
      for (j in (i + 1):K) {
        idxcount <- idxcount + 1
        b_i <- b[((i - 1) * L + 1):(i * L), , drop = FALSE]
        b_j <- b[((j - 1) * L + 1):(j * L), , drop = FALSE]
        v[idxcount, ] <- v[idxcount, ] + theta * (as.vector(b_i - b_j) - alpha[idxcount, ])
      }
    }
    alpha_vector <- as.vector(t(alpha))
    v_vector <- as.vector(t(v))
    matrix_to_invert <- (1 / K) * crossprod(bar_Gamma) + phi * R + theta * crossprod(Omega)
    diag(matrix_to_invert) <- diag(matrix_to_invert) + 1e-3
    inv_matrix <- chol2inv(chol(matrix_to_invert))
    b <- inv_matrix %*% (
      (1 / K) * t(bar_Gamma) %*% bar_Y +
        theta * t(Omega) %*% (alpha_vector - v_vector / theta)
    )
    r <- Omega %*% b - alpha_vector
    if (sqrt(sum(r^2)) < 1e-4) break
  }
  list(b = b, alpha = alpha_vector, v = v_vector, basis = basis)
}

fuse_group <- function(beta.basis, b.mat, tol = 0) {
  Z <- inprod(beta.basis, beta.basis)
  if (tol != 0) {
    tau <- tol
  } else {
    tau <- sqrt(mean(rowSums((b.mat %*% Z) * b.mat))) * 0.5
  }
  grp <- 1:nrow(b.mat)
  for (i in 1:(nrow(b.mat) - 1)) {
    for (j in (i + 1):nrow(b.mat)) {
      tmp <- b.mat[i, ] - b.mat[j, ]
      difference <- sqrt(t(tmp) %*% Z %*% tmp)
      if (difference <= tau) {
        grp[j] <- grp[i]
        b.mat[j, ] <- b.mat[i, ]
      }
    }
  }
  list(group = grp, b.hat = b.mat)
}

