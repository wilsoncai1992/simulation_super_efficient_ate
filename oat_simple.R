oat_intial_fit <- function(data_sim, is_oat = TRUE) {
  df_x <- data.frame(A = data_sim$A, W = data_sim$W)
  # weights, lambda
  hal_q <- hal9001::fit_hal(
    X = df_x,
    Y = data_sim$Y,
    family = "gaussian",
    fit_type = "glmnet",
    n_folds = 3,
    use_min = TRUE,
    yolo = FALSE,
    keep = TRUE
  )

  qn_1w <- predict(hal_q, new_data = data.frame(A = 1, W = data_sim$W))
  qn_0w <- predict(hal_q, new_data = data.frame(A = 0, W = data_sim$W))
  qn_aw <- predict(hal_q, new_data = df_x)

  if (is_oat) {
    df_w_g <- data.frame(qn_1w, qn_0w)
    hal_g <- hal9001::fit_hal(
      X = df_w_g,
      Y = data_sim$A,
      family = "binomial",
      fit_type = "glmnet",
      n_folds = 3,
      use_min = TRUE,
      yolo = FALSE,
      keep = TRUE
    )
    gn_w <- predict(hal_g, new_data = df_w_g)
  } else {
    hal_g <- hal9001::fit_hal(
      X = data_sim$W,
      Y = data_sim$A,
      family = "binomial",
      fit_type = "glmnet",
      n_folds = 3,
      use_min = TRUE,
      yolo = FALSE,
      keep = TRUE
    )
    gn_w <- predict(hal_g, new_data = data_sim$W)
  }
  return(list(
    qn_0w = qn_0w,
    qn_1w = qn_1w,
    qn_aw = qn_aw,
    gn_w = gn_w,
    hal_q = hal_q,
    hal_g = hal_g
  ))
}

truth_fit <- function(data_sim) {
  qn_1w <- data_sim$Q1(data_sim$W[, 1])
  qn_0w <- data_sim$Q0(data_sim$W[, 1])
  qn_aw <- data_sim$A * qn_1w + (1 - data_sim$A) * qn_0w
  gn_w <- data_sim$pw(data_sim$W[, 1])
  return(list(
    qn_0w = qn_0w,
    qn_1w = qn_1w,
    qn_aw = qn_aw,
    gn_w = gn_w
  ))
}

fit_tmle_ate <- function(data_sim, is_oat = TRUE, is_oracle = FALSE) {
  if (is_oracle) {
    hal_fit <- truth_fit(data_sim)
  } else {
    hal_fit <- oat_intial_fit(data_sim, is_oat = is_oat)
  }
  qn_0w <- hal_fit$qn_0w
  qn_1w <- hal_fit$qn_1w
  # qn_aw <- hal_fit$qn_aw
  gn_w <- hal_fit$gn_w
  # hal_q <- hal_fit$hal_q
  # hal_g <- hal_fit$hal_g

  targeted <- tmle::tmle(
    Y = data_sim$Y,
    A = data_sim$A,
    W = as.matrix(data_sim$W),
    Q = cbind(qn_0w, qn_1w),
    g1W = gn_w,
    family = "gaussian",
    fluctuation = "logistic",
    V = 3,
    verbose = FALSE
  )
  output <- targeted$estimates$ATE
  names(output) <- c("Psi", "var", "CI", "p_value")
  return(output)
}

oracle_ci_tmle_ate <- function(data_sim) {
  return(fit_tmle_ate(data_sim, is_oracle = TRUE))
}

compute_metric <- function(result, Psi_0, method) {
  bias <- result$Psi - Psi_0
  mse <- (result$Psi - Psi_0) ^ 2
  is_cover <- (result$CI[1] <= Psi_0) & (result$CI[2] >= Psi_0)
  output <- data.frame(
    method = method,
    bias = bias,
    mse = mse,
    is_cover = is_cover
  )
  return(output)
}
