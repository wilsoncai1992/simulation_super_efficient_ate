library(hal9001)
library(tmle)
library(dplyr)

source("./david_oat.R")
source("./david_ks.R")

N_SIMULATION <- 5e2
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

Psi_0 <- 0
# n_sim <- 1e3
# n_grid <- c(1e3)
n_grid <- c(1e2, 1e3)
# n_grid <- c(1e2, 1e3, 1e4)
df_simulation_result <- foreach(
  n_sim = n_grid,
  .combine = rbind,
  .packages = c("data.table", "drtmle", "hal9001", "SuperLearner"),
  .inorder = FALSE,
  .errorhandling = "pass",
  .verbose = TRUE
) %:%
  foreach(it2 = 1:N_SIMULATION, .combine = rbind) %dopar% {
    do_once(1, n = n_sim)
  }
head(df_simulation_result)

create_result_for_oracle_ci <- function(df_simulation_result, method = "oat") {
  df_result <- df_simulation_result[df_simulation_result$method == method, ]
  df_result <- df_result %>% group_by(n) %>% mutate(sd = sd(bias))
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

df_mc_result <- df_simulation_result %>%
  group_by(method, n) %>%
  summarize(
    bias = mean(bias),
    mse = mean(mse),
    variance = mse - bias ^ 2,
    coverage = mean(is_cover),
    sd = mean(sd),
    count = dplyr::n()
  )

save(
  df_mc_result,
  df_simulation_result,
  file = paste("df_mc_result.rda", sep = "")
)
closeCluster(cl)
# stopCluster(cl)
