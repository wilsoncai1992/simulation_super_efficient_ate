fit_onestep <- function(data_sim, is_oat = FALSE, is_oracle = FALSE, cv_variance = FALSE) {
  n <- length(data_sim$A)
  if (is_oracle) {
    hal_fit <- truth_fit(data_sim)
  } else {
    hal_fit <- oat_intial_fit(data_sim, is_oat = is_oat)
  }
  qn_0w <- hal_fit$qn_0w
  qn_1w <- hal_fit$qn_1w
  qn_aw <- hal_fit$qn_aw
  gn_w <- hal_fit$gn_w
  # hal_q <- hal_fit$hal_q
  # hal_g <- hal_fit$hal_g
  lambda_q <- hal_fit$hal_q$lambda_star
  lambda_g <- hal_fit$hal_g$lambda_star

  compute_eic <- function(A, gk, Y, Qk, Q1k, Q0k, psi) {
    HA <- A / gk - (1 - A) / (1 - gk)
    EIC <- HA * (Y - Qk) + Q1k - Q0k - psi
    return(EIC)
  }

  psi_simple <- mean(qn_1w - qn_0w)
  eic_n <- compute_eic(
    A = data_sim$A,
    gk = gn_w,
    Y = data_sim$Y,
    Qk = qn_aw,
    Q1k = qn_1w,
    Q0k = qn_0w,
    psi = psi_simple
  )
  psi_n <- psi_simple + mean(eic_n)
  var_n <- var(eic_n) / n
  se_n <- sqrt(var_n)
  CI <- psi_n + c(-1.96, 1.96) * se_n
  p_value <- 2 * pnorm(abs(psi_n / se_n), lower.tail = FALSE)

  output <- list(Psi = psi_n, var = var_n, CI = CI, p_value = p_value)
  if (cv_variance) {
    fit_prevalidated <- get_prevalidated_array(
      data_sim,
      lambda_q,
      lambda_g,
      V = 5,
      is_oat = is_oat,
      hal_fit = hal_fit
    )
    psi_simple_2 <- mean(fit_prevalidated$qn_1w - fit_prevalidated$qn_0w)
    eic_n_2 <- compute_eic(
      A = data_sim$A,
      gk = fit_prevalidated$gn_w,
      Y = data_sim$Y,
      Qk = fit_prevalidated$qn_aw,
      Q1k = fit_prevalidated$qn_1w,
      Q0k = fit_prevalidated$qn_0w,
      psi = psi_simple_2
    )
    psi_n_2 <- psi_simple_2 + mean(eic_n_2)
    var_n_2 <- var(eic_n_2) / n
    se_n_2 <- sqrt(var_n_2)
    # CI_2 <- psi_n_2 + c(-1.96, 1.96) * se_n_2
    # p_value <- 2 * pnorm(abs(psi_n_2 / se_n_2), lower.tail = FALSE)
    CI_2 <- psi_n + c(-1.96, 1.96) * se_n_2
    p_value_2 <- 2 * pnorm(abs(psi_n / se_n_2), lower.tail = FALSE)
    output <- list(Psi = psi_n, var = var_n_2, CI = CI_2, p_value = p_value_2)
  }
  return(output)
}

oracle_ci_onestep <- function(data_sim) {
  return(fit_onestep(data_sim, is_oat = FALSE, is_oracle = TRUE, cv_variance = FALSE))
}
