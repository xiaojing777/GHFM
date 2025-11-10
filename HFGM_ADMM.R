library(fda)
library(Matrix)

compute_basis_custom <- function(T, d, M) {
  create.bspline.basis(rangeval = c(0, T), nbasis = M + d, norder = d + 1)
}

compute_R_custom <- function(basis, n) {
  R1 <- eval.penalty(basis, int2Lfd(2))
  rlist <- replicate(n, R1, simplify = FALSE)
  R_block <- bdiag(rlist)
  as.matrix(R_block)
}

compute_Gamma_custom <- function(X.fd, basis) {
  n <- ncol(X.fd$coefs)
  X.coef <- t(X.fd$coefs)
  inprod_mat <- inprod(X.fd$basis, basis)
  Gamma <- X.coef %*% inprod_mat
  Gamma
}

compute_H_block <- function(Gamma) {
  n <- nrow(Gamma)
  L <- ncol(Gamma)
  H_list <- vector("list", n)
  for (i in 1:n) {
    mat <- matrix(0, nrow = 1, ncol = n * L)
    mat[1, ((i - 1) * L + 1):(i * L)] <- Gamma[i, ]
    H_list[[i]] <- mat
  }
  H <- do.call(rbind, H_list)
  H
}

compute_H_block_with_pre_grps <- function(Gamma, pre_grp) {
  K <- length(unique(pre_grp))
  n <- nrow(Gamma)
  L <- ncol(Gamma)
  H_list <- vector("list", n)
  for (k in 1:K) {
    idx.this.grp <- which(pre_grp == k)
    for (i in idx.this.grp) {
      mat <- matrix(0, nrow = 1, ncol = K * L)
      mat[1, (L * k - L + 1):(L * k)] <- Gamma[i, ]
      H_list[[i]] <- mat
    }
  }
  H <- do.call(rbind, H_list)
  H
}

compute_Omega_custom <- function(n, L) {
  cmat <- combn(n, 2)
  nrow_epsilon <- ncol(cmat)
  row_inds <- rep(seq_len(nrow_epsilon), each = 2)
  col_inds <- as.vector(rbind(cmat[1, ], cmat[2, ]))
  vals <- as.vector(rbind(rep(1, nrow_epsilon), rep(-1, nrow_epsilon)))
  eps_sparse <- sparseMatrix(i = row_inds, j = col_inds, x = vals, dims = c(nrow_epsilon, n))
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

hfgm_admm_without_pre_grp <- function(X.fd, Y, T, d, M, phi, lambda, theta, tau1, tau2, max_iter = 100, init_b = NULL) {
  Y <- matrix(Y, ncol = 1)
  n <- nrow(Y)
  basis <- compute_basis_custom(T, d, M)
  L <- basis$nbasis
  R <- compute_R_custom(basis, n)
  Gamma <- compute_Gamma_custom(X.fd, basis)
  H <- compute_H_block(Gamma)
  Omega <- compute_Omega_custom(n, L)
  matrix_to_invert <- (1 / n) * crossprod(H) + phi * R + theta * crossprod(Omega)
  diag(matrix_to_invert) <- diag(matrix_to_invert) + 1e-3
  chol_matrix <- chol(matrix_to_invert)
  inv_matrix <- chol2inv(chol_matrix)
  if (is.null(init_b)) {
    b <- inv_matrix %*% ((1 / n) * t(H) %*% Y)
  } else {
    b <- init_b
  }
  b <- as.matrix(b)
  alpha <- matrix(0, nrow = n * (n - 1) / 2, ncol = L)
  v <- matrix(0, nrow = n * (n - 1) / 2, ncol = L)
  for (m in 1:max_iter) {
    idxcount <- 0
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        idxcount <- idxcount + 1
        b_i <- b[((i - 1) * L + 1):(i * L), ]
        b_j <- b[((j - 1) * L + 1):(j * L), ]
        diff_b <- b_i - b_j
        alpha[idxcount, ] <- shrinkage(diff_b + v[idxcount, ] / theta, lambda / theta)
      }
    }
    idxcount <- 0
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        idxcount <- idxcount + 1
        b_i <- b[((i - 1) * L + 1):(i * L), ]
        b_j <- b[((j - 1) * L + 1):(j * L), ]
        v[idxcount, ] <- v[idxcount, ] + theta * (b_i - b_j - alpha[idxcount, ])
      }
    }
    alpha_vector <- as.vector(t(alpha))
    v_vector <- as.vector(t(v))
    b <- inv_matrix %*% ((1 / n) * t(H) %*% Y + theta * t(Omega) %*% (alpha_vector - v_vector / theta))
    r <- Omega %*% b - alpha_vector
    if (sqrt(sum(r^2)) < tau2) {
      cat("Converged at iteration", m, "\n")
      break
    }
  }
  list(b = b, alpha = alpha_vector, v = v_vector, basis = basis)
}

hfgm_admm_with_pre_grp <- function(X.fd, Y, T, d, M, phi, lambda, theta, tau1, tau2, max_iter = 100, init_b = NULL, pre_grp = NULL) {
  Y <- matrix(Y, ncol = 1)
  n <- nrow(Y)
  basis <- compute_basis_custom(T, d, M)
  L <- basis$nbasis
  K <- length(unique(pre_grp))
  R <- compute_R_custom(basis, K)
  Gamma <- compute_Gamma_custom(X.fd, basis)
  H <- compute_H_block_with_pre_grps(Gamma, pre_grp)
  Omega <- compute_Omega_custom(K, L)
  matrix_to_invert <- (1 / n) * crossprod(H) + phi * R + theta * crossprod(Omega)
  diag(matrix_to_invert) <- diag(matrix_to_invert) + 1e-3
  chol_matrix <- chol(matrix_to_invert)
  inv_matrix <- chol2inv(chol_matrix)
  if (is.null(init_b)) {
    b <- inv_matrix %*% ((1 / n) * t(H) %*% Y)
  } else {
    b <- init_b
  }
  b <- as.matrix(b)
  alpha <- matrix(0, nrow = K * (K - 1) / 2, ncol = L)
  v <- matrix(0, nrow = K * (K - 1) / 2, ncol = L)
  for (m in 1:max_iter) {
    idxcount <- 0
    for (i in 1:(K - 1)) {
      for (j in (i + 1):K) {
        idxcount <- idxcount + 1
        b_i <- b[((i - 1) * L + 1):(i * L), ]
        b_j <- b[((j - 1) * L + 1):(j * L), ]
        diff_b <- b_i - b_j
        alpha[idxcount, ] <- shrinkage(diff_b + v[idxcount, ] / theta, lambda / theta)
      }
    }
    idxcount <- 0
    for (i in 1:(K - 1)) {
      for (j in (i + 1):K) {
        idxcount <- idxcount + 1
        b_i <- b[((i - 1) * L + 1):(i * L), ]
        b_j <- b[((j - 1) * L + 1):(j * L), ]
        v[idxcount, ] <- v[idxcount, ] + theta * (b_i - b_j - alpha[idxcount, ])
      }
    }
    alpha_vector <- as.vector(t(alpha))
    v_vector <- as.vector(t(v))
    b <- inv_matrix %*% ((1 / n) * t(H) %*% Y + theta * t(Omega) %*% (alpha_vector - v_vector / theta))
    r <- Omega %*% b - alpha_vector
    if (sqrt(sum(r^2)) < tau2) {
      cat("Converged at iteration", m, "\n")
      break
    }
  }
  list(b = b, alpha = alpha_vector, v = v_vector, basis = basis)
}

fuse_group <- function(b, basis, tol = NULL) {
  L <- basis$nbasis
  n <- length(b) / L
  B_mat <- matrix(b, nrow = n, byrow = TRUE)
  Z <- inprod(basis, basis)
  if (is.null(tol)) {
    tol <- sqrt(mean(rowSums((B_mat %*% Z) * B_mat))) * 0.5
  }
  cat("Using tol =", tol, "\n")
  grp <- 1:n
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- B_mat[i, ] - B_mat[j, ]
      dist <- sqrt(t(tmp) %*% Z %*% tmp)
      if (dist <= tol) {
        grp[j] <- grp[i]
        B_mat[j, ] <- B_mat[i, ]
      }
    }
  }
  list(group = grp, b.hat = B_mat)
}

