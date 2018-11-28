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

get_prevalidated_array <- function(
  data_sim,
  lambda_q,
  lambda_g,
  V,
  is_oat = FALSE,
  hal_fit = NULL
) {
  n <- length(data_sim$A)
  generate_fold_id <- function(x, n) split(x, cut(
    sample(seq_along(x), size = length(x), replace = FALSE),
    n,
    labels = FALSE
  ))
  fold_idx <- generate_fold_id(seq_len(n), V)
  qn_1ws <- c()
  qn_0ws <- c()
  qn_aws <- c()
  gn_ws <- c()
  for (v in 1:V) {
    fold_id_valid <- fold_idx[[v]]
    fold_id_train <- unlist(fold_idx[setdiff(1:V, v)])
    names(fold_id_train) <- NULL

    data_valid <- subset_data_object(data_sim, fold_id_valid)
    data_train <- subset_data_object(data_sim, fold_id_train)
    df_x <- data.frame(A = data_train$A, W = data_train$W)
    hal_q <- hal9001::fit_hal_single_lambda(
      X = df_x,
      Y = data_train$Y,
      family = "gaussian",
      fit_type = "glmnet",
      n_folds = 3,
      use_min = TRUE,
      yolo = FALSE,
      lambda = lambda_q
    )
    qn_1w <- predict(hal_q, new_data = data.frame(A = 1, W = data_valid$W))
    qn_0w <- predict(hal_q, new_data = data.frame(A = 0, W = data_valid$W))
    qn_aw <- predict(hal_q, new_data = data.frame(A = data_valid$A, W = data_valid$W))
    qn_1ws <- c(qn_1ws, qn_1w)
    qn_0ws <- c(qn_0ws, qn_0w)
    qn_aws <- c(qn_aws, qn_aw)

    # fit g
    if (is_oat) {
      # use training data from the whole fit of the training data
      df_w_g <- data.frame(hal_fit$qn_1w, hal_fit$qn_0w)
      Y_g <- data_sim$A

      # use training data from just the fit on the training fold
      # qn_1w_train <- predict(hal_q, new_data = data.frame(A = 1, W = data_train$W))
      # qn_0w_train <- predict(hal_q, new_data = data.frame(A = 0, W = data_train$W))
      # qn_aw_train <- predict(hal_q, new_data = data.frame(A = data_train$A, W = data_train$W))
      # df_w_g <- data.frame(qn_1w_train, qn_0w_train)
      # Y_g <- data_train$A

      hal_g <- hal9001::fit_hal_single_lambda(
        X = df_w_g,
        Y = Y_g,
        family = "binomial",
        fit_type = "glmnet",
        n_folds = 3,
        use_min = TRUE,
        yolo = FALSE,
        lambda = lambda_g
      )
      # predict on the validation fold data
      gn_w <- predict(hal_g, new_data = data.frame(qn_1w, qn_0w))
    } else {
      hal_g <- hal9001::fit_hal_single_lambda(
        X = data_train$W,
        Y = data_train$A,
        family = "binomial",
        fit_type = "glmnet",
        n_folds = 3,
        use_min = TRUE,
        yolo = FALSE,
        lambda = lambda_g
      )
      gn_w <- predict(hal_g, new_data = data_valid$W)
    }
    gn_ws <- c(gn_ws, gn_w)
  }
  return(list(
    qn_1w = qn_1ws,
    qn_0w = qn_0ws,
    qn_aw = qn_aws,
    gn_w = gn_ws
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

fit_tmle_ate <- function(
  data_sim,
  is_oat = FALSE,
  is_oracle = FALSE,
  cv_variance = FALSE
) {
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
  lambda_q <- hal_fit$hal_q$lambda_star
  lambda_g <- hal_fit$hal_g$lambda_star

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

  if (cv_variance) {
    fit_prevalidated <- get_prevalidated_array(
      data_sim,
      lambda_q,
      lambda_g,
      V = 5,
      is_oat = is_oat,
      hal_fit = hal_fit
    )
    tmle_for_cv_variance <- tmle::tmle(
      Y = data_sim$Y,
      A = data_sim$A,
      W = as.matrix(data_sim$W),
      Q = cbind(fit_prevalidated$qn_0w, fit_prevalidated$qn_1w),
      g1W = fit_prevalidated$gn_w,
      family = "gaussian",
      fluctuation = "logistic",
      V = 3,
      verbose = FALSE
    )
    output_cv <- tmle_for_cv_variance$estimates$ATE
    names(output_cv) <- c("Psi", "var", "CI", "p_value")
    output$var <- output_cv$var
    output$CI <- output$Psi + c(-1.96, 1.96) * sqrt(output_cv$var)
    output$p_value <- 2 * pnorm(
      abs(output$Psi / sqrt(output_cv$var)),
      lower.tail = FALSE
    )
  }
  return(output)
}

oracle_ci_tmle_ate <- function(data_sim) {
  return(fit_tmle_ate(data_sim, is_oat = TRUE, is_oracle = TRUE, cv_variance = FALSE))
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


subset_data_object <- function(data, index) {
  data$Y <- data$Y[index]
  data$A <- data$A[index]
  data$W <- data$W[index, ]
  data$X_matrix <- data$X_matrix[index, ]
  data$all_df <- data$all_df[index, ]
  return(data)
}
