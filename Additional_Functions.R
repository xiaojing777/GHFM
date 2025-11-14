compute_ISE <- function(bhat_mat, btrue_mat, G) {
  diff_mat <- bhat_mat - btrue_mat
  err_sq   <- sum((diff_mat %*% G) * diff_mat)
  true_sq  <- sum((btrue_mat %*% G) * btrue_mat)
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
