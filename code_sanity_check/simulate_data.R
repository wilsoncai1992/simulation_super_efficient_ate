simulate_data <- function(n_sim, a1, a2, b1) {
  thresholding <- function(x, min, max) pmin(pmax(x, min), max)

  W <- truncnorm::rtruncnorm(n = n_sim, a = -10, b = 10, mean = 0, sd = 4)
  A <- rbinom(
    n_sim,
    size = 1,
    prob = thresholding(.3 + 0.1 * W * sin(a2 * W), 0.3, 0.7) +
      rnorm(n_sim, mean = 0, sd = 0.05)
  )

  Y <- b1 * sin(W * a1) + A + rnorm(n_sim, 0, 1)
  Q1 <- function(w) return(b1 * sin(w * a1) + 1)
  Q0 <- function(w) return(b1 * sin(w * a1))
  pw <- function(w) return(thresholding(.3 + 0.1 * w * sin(a2 * w), 0.3, 0.7))

  X_matrix_0 <- data.frame(A, W)
  all_df <- data.frame(Y, A, W)

  # append one value of Z
  Z <- rep(1, n_sim)
  X_matrix <- cbind(Z, X_matrix_0)
  all_df <- cbind(Z, all_df)
  return(list(
    Y = Y,
    A = A,
    W = data.frame(W),
    X_matrix = X_matrix,
    all_df = all_df,
    Q1 = Q1,
    Q0 = Q0,
    pw = pw
  ))
}
