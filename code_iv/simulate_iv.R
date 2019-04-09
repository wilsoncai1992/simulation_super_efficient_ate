# DGP:
# W_1,...W_{n_covar} ~ U_{n_covar}(-1.5, 1.5)
# W_{n_covar + 1} is the instrument with Bernoulli(1/2) dist.
# \bar{G}_0(W) = expit(iv_beta/2 - iv_beta * W_{n_covar + 1} + sum_{j=1}^{n_covar} \beta_j W_j ) ,
# where \beta_j = 2^{1 - j}, j = 1,...,n_covar
# For practical purposes, we keep drawing A, until we get at least some observations
# with instrument = 0 and A = 0/1, and instrument = 1 and A = 0/1
# \bar{Q}_0(A,W) = A - sum_{j=1}^{n_covar} \beta_j W_j , where \beta_j are as above
# and notice that instrument has no effect on \bar{Q}
# Y = \bar{Q}_0(A,W) + N(0,1)

simulate_iv <- function(n, n_covar = 8, iv_beta = 12) {
  W <- matrix(runif((n_covar - 1) * n, -1.5, 1.5), ncol = n_covar - 1)
  # the instrument
  W <- cbind(W, rbinom(n, 1, 1 / 2))

  # find the max and min propensity score (to report in text)
  # W_max <- matrix(1.5, nrow = 1, ncol = 7)
  # W_max <- cbind(W_max, 0)

  # W_min <- matrix(-1.5, nrow = 1, ncol = 7)
  # W_min <- cbind(W_min, 1)

  # lg0 <- apply(W_min, 1, function(w){
  #     iv_beta/2 - iv_beta*w[length(w)] + sum(2^(-(0:(n_covar-2))) * w[-length(w)])
  # })
  # plogis(lg0)

  # lg0 <- apply(W_max, 1, function(w){
  #     iv_beta/2 - iv_beta*w[length(w)] + sum(2^(-(0:(n_covar-2))) * w[-length(w)])
  # })
  # plogis(lg0)

  logit_G0 <- apply(W, 1, function(w) {
    iv_beta / 2 - iv_beta * w[length(w)] + sum(2^(-(0:(n_covar - 2))) * w[-length(w)])
  })

  G0 <- plogis(logit_G0)

  # make sure at least some overlap
  cross_tab <- rep(0, 4)
  while (any(cross_tab == 0)) {
    A <- rbinom(n, 1, G0)
    cross_tab <- as.numeric(table(A, W[, ncol(W)]))
  }

  QAW <- apply(W, 1, function(w) {
    -sum(2^(-(0:(n_covar - 2))) * w[-length(w)])
  })
  QAW <- QAW + A
  Q1W <- apply(W, 1, function(w) {
    1 - sum(2^(-(0:(n_covar - 2))) * w[-length(w)])
  })
  Q0W <- apply(W, 1, function(w) {
    0 - sum(2^(-(0:(n_covar - 2))) * w[-length(w)])
  })

  Y <- QAW + rnorm(n)
  W_df <- data.frame(W)
  colnames(W_df) <- paste0("W", 1:n_covar)

  Q1W <- function(W, n_covar = n_covar) {
    apply(W, 1, function(w) {
      1 - sum(2^(-(0:(n_covar - 2))) * w[-length(w)])
    })
  }

  Q0W <- function(W, n_covar = n_covar) {
    apply(W, 1, function(w) {
      0 - sum(2^(-(0:(n_covar - 2))) * w[-length(w)])
    })
  }

  G <- function(W, iv_beta = iv_beta, n_covar = n_covar) {
    logit_G0 <- apply(W, 1, function(w) {
      iv_beta / 2 - iv_beta * w[length(w)] + sum(2^(-(0:(n_covar - 2))) * w[-length(w)])
    })
    return(plogis(logit_G0))
  }
  out <- list(
    W = W_df,
    A = A,
    Y = Y,
    Q1W = Q1W,
    Q0W = Q0W,
    G = G
  )
  return(out)
}
