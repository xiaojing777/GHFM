library(fda)

fit_full_vs_reduced_gauss <- function(Z, Y, groups) {
  n <- nrow(Z); L <- ncol(Z)
  df_red <- data.frame(Y = Y, Z)
  names(df_red) <- c("Y", paste0("Z", 1:L))
  form_red <- as.formula(paste("Y ~", paste(names(df_red)[-1], collapse = " + ")))
  fit_red <- lm(form_red, data = df_red)
  RSS_red <- sum(resid(fit_red)^2)
  K <- length(unique(groups))
  df_full <- data.frame(Y = Y)
  for (k in 1:K) {
    mask <- (groups == k)
    for (l in 1:L) {
      df_full[[paste0("G", k, "_Z", l)]] <- ifelse(mask, Z[, l], 0)
    }
  }
  form_full <- as.formula(paste("Y ~", paste(names(df_full)[-1], collapse = " + ")))
  fit_full <- lm(form_full, data = df_full)
  RSS_full <- sum(resid(fit_full)^2)
  df_num <- 2 * L
  df_full_err <- n - (1 + 3 * L)
  F_stat <- ((RSS_red - RSS_full) / df_num) / (RSS_full / df_full_err)
  p_val <- 1 - pf(F_stat, df_num, df_full_err)
  list(F_stat = F_stat,
       p_value = p_val,
       RSS_full = RSS_full,
       RSS_reduced = RSS_red)
}

fit_full_vs_reduced_ber <- function(Z, Y, groups) {
  n <- nrow(Z); L <- ncol(Z)
  df_red <- data.frame(Y = Y, Z)
  names(df_red) <- c("Y", paste0("Z", 1:L))
  form_red <- as.formula(paste("Y ~", paste(names(df_red)[-1], collapse = " + ")))
  fit_red <- glm(form_red, data = df_red, family = binomial(link = "logit"))
  D_red <- deviance(fit_red)
  df_red_p <- df.residual(fit_red)
  K <- length(unique(groups))
  df_full <- data.frame(Y = Y)
  for (k in 1:K) {
    mask <- (groups == k)
    for (l in 1:L) {
      df_full[[paste0("G", k, "_Z", l)]] <- ifelse(mask, Z[, l], 0)
    }
  }
  form_full <- as.formula(paste("Y ~", paste(names(df_full)[-1], collapse = " + ")))
  fit_full <- glm(form_full, data = df_full, family = binomial(link = "logit"))
  D_full <- deviance(fit_full)
  df_full_p <- df.residual(fit_full)
  D_diff <- D_red - D_full
  df_diff <- df_red_p - df_full_p
  p_value <- 1 - pchisq(D_diff, df_diff)
  list(deviance_reduced = D_red,
       deviance_full = D_full,
       chi_sq = D_diff,
       df = df_diff,
       p_value = p_value)
}

