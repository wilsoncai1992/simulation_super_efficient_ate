library(dplyr)
library(hal9001)
library(glmnet)

fit_nuisance <- function(W, A, Y){
  Q_fit <- cv.glmnet(
    x = cbind(A, W), y = Y, nfolds = 5, family = "binomial", standardize = FALSE
  )
  Q1W <- predict(Q_fit, newx = cbind(1, W), s = "lambda.min", type = "response")[, 1]
  Q0W <- predict(Q_fit, newx = cbind(0, W), s = "lambda.min", type = "response")[, 1]
  QAW <- predict(Q_fit, newx = cbind(A, W), s = "lambda.min", type = "response")[, 1]

  G_fit <- cv.glmnet(
    x = W, y = A, nfolds = 5, family = "binomial", standardize = FALSE
  )
  G1W <- predict(G_fit, newx = W, s = "lambda.min", type = "response")[, 1]
  return(data.frame(Q1W, Q0W, QAW, G1W))
}
compute_eic <- function(Q1W, Q0W, G1W, Y, A, psi_n) {
  return((2 * A - 1) / G1W * (Y - Q1W + Q0W) + (mean(Q1W - Q0W) - psi_n))
}

fit_one_A <- function(W, A, Y, delta = 0.1) {
  nuisance_fits <- fit_nuisance(W = W, A = A, Y = Y)
  propens_not_extreme <- mean(
    (nuisance_fits$G1W >= delta) & (nuisance_fits$G1W <= (1 - delta))
  )

  feature_reduced <- as.matrix(
    data.frame(Q1W = nuisance_fits$Q1W, Q0W = nuisance_fits$Q0W)
  )
  if (isTRUE(all.equal(feature_reduced[,1], feature_reduced[,2]))) {
    feature_reduced <- matrix(feature_reduced[, 1])
  }
  if (min(feature_reduced) == max(feature_reduced)) {
    # if the reduced Q fits are horizonal; just fit a global mean G fit
    nuisance_fits$G1W_reduced <- mean(A)
  } else {
    G_fit_reduced <- hal9001::fit_hal(
      X = feature_reduced,
      Y = A,
      fit_type = 'glmnet',
      n_folds = 5,
      use_min = TRUE,
      family = 'binomial',
      return_lasso = TRUE,
      cv_select = TRUE,
      yolo = FALSE
    )
    nuisance_fits$G1W_reduced <- predict(
      G_fit_reduced, new_data = feature_reduced, type = "response"
    )
  }

  tmle_regular <- tmle::tmle(
    Y = Y,
    A = A,
    W = W,
    Q = as.matrix(data.frame(nuisance_fits$Q0W, nuisance_fits$Q1W)),
    g1W = nuisance_fits$G1W,
    gbound = 1e-2,
    family = 'binomial'
  )
  tmle_reduced <- tmle::tmle(
    Y = Y,
    A = A,
    W = W,
    Q = as.matrix(data.frame(nuisance_fits$Q0W, nuisance_fits$Q1W)),
    g1W = nuisance_fits$G1W_reduced,
    gbound = 1e-2,
    family = 'binomial'
  )
  df_tmle_regular <- data.frame(
    psi = tmle_regular$estimates$ATE$psi,
    std_err = sqrt(tmle_regular$estimates$ATE$var.psi)
    # ci_lower = tmle_regular$estimates$ATE$CI[1],
    # ci_upper = tmle_regular$estimates$ATE$CI[2]
  )
  df_tmle_reduced <- data.frame(
    psi = tmle_reduced$estimates$ATE$psi,
    std_err = sqrt(tmle_reduced$estimates$ATE$var.psi)
    # ci_lower = tmle_reduced$estimates$ATE$CI[1],
    # ci_upper = tmle_reduced$estimates$ATE$CI[2]
  )

  psi_n <- mean(nuisance_fits$Q1W - nuisance_fits$Q0W)
  n <- nrow(nuisance_fits)
  eic_fit <- compute_eic(
    Q1W = nuisance_fits$Q1W,
    Q0W = nuisance_fits$Q0W,
    G1W = nuisance_fits$G1W,
    Y = Y,
    A = A,
    psi_n = psi_n
  )
  df_onestep_regular <- data.frame(
    psi = psi_n + mean(eic_fit),
    std_err = sd(eic_fit) / sqrt(n)
  )
  df_plugin <- data.frame(
    psi = psi_n,
    std_err = sd(eic_fit) / sqrt(n)
  )

  eic_fit_reduced <- compute_eic(
    Q1W = nuisance_fits$Q1W,
    Q0W = nuisance_fits$Q0W,
    G1W = nuisance_fits$G1W_reduced,
    Y = Y,
    A = A,
    psi_n = psi_n
  )
  df_onestep_reduced <- data.frame(
    psi = psi_n + mean(eic_fit_reduced),
    std_err = sd(eic_fit_reduced) / sqrt(n)
  )

  df_tmle_regular$method <- "tmle_regular"
  df_tmle_reduced$method <- "tmle_reduced"
  df_onestep_regular$method <- "onestep_regular"
  df_onestep_reduced$method <- "onestep_reduced"
  df_plugin$method <- "plugin"
  df_result <- rbind(
    df_tmle_regular, df_tmle_reduced, df_onestep_regular, df_onestep_reduced, df_plugin
  )
  df_result <- df_result %>%
    mutate(ci_lower = psi - 1.96 * std_err, ci_upper = psi + 1.96 * std_err)
  df_result$A <- A_name
  df_result$positivity_score <- propens_not_extreme
  return(df_result)
}

