library(hal9001)
library(tmle)
library(dplyr)

source("./david_oat_sl.R")
source("./david_realdata.R")

N_SIMULATION <- 2e2
library(foreach)
library(Rmpi)
library(doMPI)
cl = startMPIcluster()
registerDoMPI(cl)
clusterSize(cl) # just to check

Psi_0 <- 0
df_simulation_result <- foreach(
  i_data = 2:2,
  .combine = c,
  .packages = c("data.table", "drtmle", "hal9001", "SuperLearner"),
  .inorder = FALSE,
  .errorhandling = "pass",
  .verbose = TRUE
) %:%
  foreach(it2 = 1:N_SIMULATION, .combine = c, .errorhandling = "remove") %dopar% {
    if (i_data == 1) {load("./fev_data.RData"); dat <- fev_data; dat_name <- "fev"}
    if (i_data == 2) {load("./cebu_data.RData"); dat <- cebu_data; dat_name <- "cebu"}
    if (i_data == 3) {load("./wine_data.RData"); dat <- wine_data; dat_name <- "wine"}
    do_once(dat, dat_name, permute_y = TRUE)
  }
head(df_simulation_result)

df_result <- df_simulation_result[names(df_simulation_result) == "df_result"]
df_result <- do.call(rbind, df_result)
Ys <- df_simulation_result[names(df_simulation_result) == "Y"]
Ys <- do.call(rbind, Ys)

save(
  df_result,
  Ys,
  file = paste("debug.rda", sep = "")
)

closeCluster(cl)
mpi.quit()
# stopCluster(cl)
