library(hal9001)
library(tmle)
library(dplyr)

source("./simulate_data.R")
source("./oat_simple.R")
source("./onestep_ate.R")

do_once <- function(n_sim) {
  data_sim <- simulate_data(n_sim = n_sim, a1 = a1, a2 = a2, b1 = b1)
  Psi_0 <- 1

  result_oat <- fit_tmle_ate(data_sim, is_oat = TRUE)
  result_onestep <- fit_onestep(data_sim, is_oat = TRUE)
  result_tmle <- fit_tmle_ate(data_sim, is_oat = FALSE)
  result_oat_oracle_ci <- oracle_ci_tmle_ate(data_sim)
  result_onestep_oracle_ci <- oracle_ci_onestep(data_sim)

  df_simulation_result <- list()
  i <- 1
  df_simulation_result[[i]] <- compute_metric(
    result_oat, Psi_0, "oat"
  ); i <- i + 1
  df_simulation_result[[i]] <- compute_metric(
    result_oat_oracle_ci, Psi_0, "oat_oracle_ci"
  ); i <- i + 1
  df_simulation_result[[i]] <- compute_metric(
    result_onestep, Psi_0, "onestep"
  ); i <- i + 1
  df_simulation_result[[i]] <- compute_metric(
    result_onestep_oracle_ci, Psi_0, "onestep_oracle_ci"
  ); i <- i + 1
  df_simulation_result[[i]] <- compute_metric(
    result_tmle, Psi_0, "tmle"
  ); i <- i + 1

  df_simulation_result <- do.call(rbind, df_simulation_result)
  df_simulation_result$n <- n_sim
  return(df_simulation_result)
}

N_SIMULATION <- 1e3
# N_SIMULATION <- 8
library(foreach)
library(Rmpi)
library(doMPI)
cl = startMPIcluster()
registerDoMPI(cl)
clusterSize(cl) # just to check

# library(doSNOW)
# library(tcltk)
# nw <- parallel:::detectCores()  # number of workers
# cl <- makeSOCKcluster(nw)
# registerDoSNOW(cl)

# n_sim <- 1e3
a1 <- .5
b1 <- 3
a2 <- .1
n_grid <- c(1e2, 1e3)
# n_grid <- c(1e2, 1e3, 1e4)
df_simulation_result <- foreach(
  n_sim = n_grid,
  .combine = rbind,
  .packages = c("R6", "tmle", "hal9001"),
  .inorder = FALSE,
  .errorhandling = "pass",
  .verbose = TRUE
) %:%
  foreach(it2 = 1:N_SIMULATION, .combine = rbind) %dopar% {
    do_once(n_sim = n_sim)
  }
head(df_simulation_result)
df_mc_result <- df_simulation_result %>%
  group_by(method, n) %>%
  summarize(
    bias = mean(bias),
    mse = mean(mse),
    variance = mse - bias ^ 2,
    coverage = mean(is_cover),
    count = dplyr::n()
  )

save(
  df_mc_result,
  df_simulation_result,
  file = paste("df_mc_result.rda", sep = "")
)
closeCluster(cl)
# stopCluster(cl)

# # prevalidated prediction
# hal_q$hal_lasso$fit.preval
# lambda_min_index <- which(hal_q$hal_lasso$lambda == hal_q$hal_lasso$lambda.min)
# yhat_validation <- hal_q$hal_lasso$fit.preval[, lambda_min_index]
