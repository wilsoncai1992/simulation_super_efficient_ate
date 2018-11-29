library(hal9001)
library(tmle)
library(dplyr)

source("./david_oat.R")
source("./simulate_data.R")
source("./david_testsim.R")


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
Psi_0 <- 1
# n_grid <- c(1e2)
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
  sd_oracle <- sd(df_result$bias)
  df_result$is_cover <- abs(df_result$bias) <= 1.96 * sd_oracle
  df_result$method <- paste(method, "_oracle_ci", sep = "")
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
    count = dplyr::n()
  )

save(
  df_mc_result,
  df_simulation_result,
  file = paste("df_mc_result.rda", sep = "")
)
closeCluster(cl)
# stopCluster(cl)
