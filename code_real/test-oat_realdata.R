library(hal9001)
library(tmle)
library(dplyr)

source("./david_oat_sl.R")
source("./david_realdata.R")

N_SIMULATION <- 1e3
# N_SIMULATION <- 2
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
df_simulation_result <- foreach(
  i_data = 1:3,
  .combine = rbind,
  .packages = c("data.table", "drtmle", "hal9001", "SuperLearner"),
  .inorder = FALSE,
  .errorhandling = "pass",
  .verbose = TRUE
) %:%
  foreach(it2 = 1:N_SIMULATION, .combine = rbind, .errorhandling = "remove") %dopar% {
    if (i_data == 1) {load("./fev_data.RData"); dat <- fev_data; dat_name <- "fev"}
    if (i_data == 2) {load("./cebu_data.RData"); dat <- cebu_data; dat_name <- "cebu"}
    if (i_data == 3) {load("./wine_data.RData"); dat <- wine_data; dat_name <- "wine"}
    do_once(dat, dat_name, permute_y = TRUE)
  }
head(df_simulation_result)

create_result_for_oracle_ci <- function(df_simulation_result, method = "oat") {
  df_result <- df_simulation_result[df_simulation_result$method == method, ]
  df_result <- df_result %>% group_by(dat_name) %>% mutate(sd = sd(bias))
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
  group_by(method, dat_name) %>%
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

df_no_permute <- list()
for (i_data in 1:3) {
  if (i_data == 1) {load("./fev_data.RData"); dat <- fev_data; dat_name <- "fev"}
  if (i_data == 2) {load("./cebu_data.RData"); dat <- cebu_data; dat_name <- "cebu"}
  if (i_data == 3) {load("./wine_data.RData"); dat <- wine_data; dat_name <- "wine"}
  df_no_permute <- c(
    df_no_permute, list(do_once(dat, dat_name, permute_y = FALSE))
  )
}
df_no_permute <- do.call(rbind, df_no_permute)
save(
  df_no_permute,
  file = paste("df_no_permute.rda", sep = "")
)

closeCluster(cl)
mpi.quit()
# stopCluster(cl)
