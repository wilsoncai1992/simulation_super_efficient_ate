library(hal9001)
library(tmle)
library(dplyr)
library(drtmle)

source("./oat_hal_ols.R")
source("./simulate_iv.R")
source("./do_once_iv.R")

N_SIMULATION <- 1e3
# N_SIMULATION <- 4
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

Psi_0 <- 1
# n_grid <- c(1e2)
n_grid <- c(1e2, 5e2, 1e3)

iv_beta_grid <- c(0, 3, 6)
df_simulation_result <- foreach(
  n_sim = n_grid,
  .combine = rbind,
  .packages = c("data.table", "drtmle", "hal9001", "SuperLearner"),
  .inorder = FALSE,
  .errorhandling = "pass",
  .verbose = TRUE
) %:%
  foreach(iv_beta = iv_beta_grid, .combine = rbind) %:%
    foreach(it2 = 1:N_SIMULATION, .combine = rbind) %dopar% {
      do_once(
      #   1, n = n_sim, n_covar = 8, iv_beta = iv_beta, full_adaptive_cv = FALSE
        1, n = n_sim, n_covar = 8, iv_beta = iv_beta, full_adaptive_cv = TRUE
      )
    }
head(df_simulation_result)

create_result_for_oracle_ci <- function(df_simulation_result, method = "oat") {
  df_result <- df_simulation_result[df_simulation_result$method == method, ]
  df_result <- df_result %>% group_by(n, iv_beta, n_covar) %>% mutate(sd = sd(bias))
  df_result$is_cover <- abs(df_result$bias) <= 1.96 * df_result$sd
  df_result$method <- paste(method, "_oracle_ci", sep = "")
  df_result <- df_result %>% ungroup()
  return(df_result)
}
df_simulation_result <- rbind(
  df_simulation_result,
  create_result_for_oracle_ci(df_simulation_result, "oat")
)
df_simulation_result <- rbind(
  df_simulation_result,
  create_result_for_oracle_ci(df_simulation_result, "onestep")
)
df_simulation_result <- rbind(
  df_simulation_result,
  create_result_for_oracle_ci(df_simulation_result, "onestep_oat")
)
df_simulation_result <- rbind(
  df_simulation_result,
  create_result_for_oracle_ci(df_simulation_result, "tmle")
)

df_mc_result <- df_simulation_result %>%
  group_by(method, n, n_covar, iv_beta) %>%
  summarize(
    bias = mean(bias),
    mse = mean(mse),
    coverage = mean(is_cover),
    sd = mean(sd),
    count = dplyr::n()
  ) %>%
  mutate(variance = mse - bias ^ 2)


save(
  df_mc_result,
  df_simulation_result,
  file = paste("df_mc_result.rda", sep = "")
)
closeCluster(cl)
mpi.quit()
# stopCluster(cl)