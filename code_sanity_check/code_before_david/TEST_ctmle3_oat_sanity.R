library(origami)
library(sl3)
library(hal9001)
library(tmle3)
library(ctmle3)

source('./simulate_data.R')
n_sim <- 1e3
a1 <- .5
b1 <- 3
a2 <- .1
data_sim <- simulate_data(n_sim = n_sim, a1 = a1, a2 = a2, b1 = b1)
df <- data_sim$all_df
node_list <- list(
  W = c("W"),
  A = "A",
  Y = "Y"
)

tmle_spec <- tmle_oat_TSM_all()
tmle_task <- tmle_spec$make_tmle_task(df, node_list)

hal_Q <- sl3::Lrnr_hal9001$new(
  fit_type = "glmnet",
  n_folds = 3,
  use_min = TRUE
)

learner_list <- list(Y = hal_Q, A = hal_Q)
oatmle_fit <- tmle3::tmle3(tmle_spec, df, node_list, learner_list = learner_list)
# extract results
tmle3_psi <- oatmle_fit$summary$tmle_est
tmle3_se <- oatmle_fit$summary$se
tmle3_epsilon <- oatmle_fit$updater$epsilons[[1]]$Y
tmle3_psi
tmle3_se
tmle3_epsilon

get_ATE <- function(tmle3_psi, tmle3_se, n) {
  Psi_n <- diff(tmle3_psi)
  Psi_n_se <- sqrt(sum(tmle3_se^2))
  Psi_n_ci <- Psi_n + c(-1.96, 1.96) * Psi_n_se / sqrt(n)
  return(list(
    Psi_n = Psi_n, 
    Psi_n_se = Psi_n_se,
    Psi_n_ci = Psi_n_ci
  ))  
}
result_ATE <- get_ATE(tmle3_psi, tmle3_se, n_sim)
Psi_n <- result_ATE$Psi_n
Psi_n_se <- result_ATE$Psi_n_se
Psi_n_ci <- result_ATE$Psi_n_ci
