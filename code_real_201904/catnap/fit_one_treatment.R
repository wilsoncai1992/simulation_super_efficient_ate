library(dplyr)
library(hal9001)
library(glmnet)
library(future)
library(listenv)
load("./all_data.RData")
A_names <- names(all_data)

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

fit_one_A <- function(W, A, Y) {
  nuisance_fits <- fit_nuisance(W = W, A = A, Y = Y)

  feature_reduced <- as.matrix(
    data.frame(Q1W = nuisance_fits$Q1W, Q0W = nuisance_fits$Q0W)
  )
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

  tmle_regular <- tmle::tmle(
    Y = Y,
    A = A,
    W = W,
    Q = as.matrix(data.frame(nuisance_fits$Q0W, nuisance_fits$Q1W)),
    g1W = nuisance_fits$G1W,
    gbound = 1e-4,
    family = 'binomial'
  )
  tmle_reduced <- tmle::tmle(
    Y = Y,
    A = A,
    W = W,
    Q = as.matrix(data.frame(nuisance_fits$Q0W, nuisance_fits$Q1W)),
    g1W = nuisance_fits$G1W_reduced,
    gbound = 1e-4,
    family = 'binomial'
  )
  df_tmle_regular <- data.frame(
    psi = tmle_regular$estimates$ATE$psi,
    std_err = sqrt(tmle_regular$estimates$ATE$var.psi),
    ci_lower = tmle_regular$estimates$ATE$CI[1],
    ci_upper = tmle_regular$estimates$ATE$CI[2]
  )
  df_tmle_reduced <- data.frame(
    psi = tmle_reduced$estimates$ATE$psi,
    std_err = sqrt(tmle_reduced$estimates$ATE$var.psi),
    ci_lower = tmle_reduced$estimates$ATE$CI[1],
    ci_upper = tmle_reduced$estimates$ATE$CI[2]
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
  df_onestep_regular$ci_lower <- df_onestep_regular$psi - 1.96 * df_onestep_regular$std_err
  df_onestep_regular$ci_upper <- df_onestep_regular$psi + 1.96 * df_onestep_regular$std_err

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
  df_onestep_reduced$ci_lower <- df_onestep_reduced$psi - 1.96 * df_onestep_reduced$std_err
  df_onestep_reduced$ci_upper <- df_onestep_reduced$psi + 1.96 * df_onestep_reduced$std_err

  df_tmle_regular$method <- "tmle_regular"
  df_tmle_reduced$method <- "tmle_reduced"
  df_onestep_regular$method <- "onestep_regular"
  df_onestep_reduced$method <- "onestep_reduced"
  df_result <- rbind(
    df_tmle_regular, df_tmle_reduced, df_onestep_regular, df_onestep_reduced
  )
  df_result$A <- A_name
  return(df_result)
}

plan(multisession, workers = 8)
# plan(sequential)
df_results <- listenv()
j <- 1
for (A_name in A_names) {
  df_results[[j]] %<-% {
    data_train <- all_data[[A_name]]
    # all W are integer type
    W <- as.matrix(data_train$W)
    A <- data_train$A
    is_not_missing <- complete.cases(W)
    W <- W[is_not_missing, ]
    A <- A[is_not_missing]
    Y <- data_train$Y[is_not_missing]
    fit_one_A(W = W, A = A, Y = Y)
  }
  j <- j + 1
}

df_results <- do.call(rbind, as.list(df_results))
table(df_results$A)

save(df_results, file = "df_results.rda")
