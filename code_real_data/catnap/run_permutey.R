load("./all_data.RData")
source("./fit_one_treatment.R")
A_names <- names(all_data)

n_mcmc <- 1e2
# n_mcmc <- 2

library(foreach)
# library(Rmpi)
# library(doMPI)
# cl = startMPIcluster()
# registerDoMPI(cl)
# clusterSize(cl) # just to check

library(doSNOW)
library(tcltk)
nw <- parallel:::detectCores() # number of workers
cl <- makeSOCKcluster(nw)
registerDoSNOW(cl)

df_results <- foreach(
  A_name = A_names,
  # A_name = A_names[1],
  .combine = rbind,
  .packages = c("R6", "hal9001", "glmnet", "dplyr"),
  .inorder = FALSE,
  .errorhandling = "remove",
  .verbose = TRUE
) %:%
  foreach(
    i = 1:n_mcmc, .combine = rbind,
    .errorhandling = "remove"
  ) %dopar% {
    data_train <- all_data[[A_name]]
    # all W are integer type
    W <- as.matrix(data_train$W)
    A <- data_train$A
    is_not_missing <- complete.cases(W)
    W <- W[is_not_missing, ]
    A <- A[is_not_missing]
    Y <- data_train$Y[is_not_missing]
    permute_y <-  TRUE
    if (permute_y) Y <- sample(Y)
    result <- fit_one_A(W = W, A = A, Y = Y)
    result$id_sim <- i
    return(result)
}


library(dplyr)
psi_0 <- 0
df_summary <- df_results %>%
  group_by(A, method) %>%
  mutate(
    bias = psi - psi_0,
    mse = (psi - psi_0) ^ 2,
    is_cover = (ci_upper >= psi_0) & (ci_lower <= psi_0)
  ) %>%
  summarise(
    bias = mean(bias),
    mse = mean(mse),
    coverage = mean(is_cover),
    cnt = length(unique(id_sim))
  ) %>%
  mutate(variance = mse - bias ^ 2)
save(df_results, df_summary, file = "permute_y.rda")

# closeCluster(cl)
# mpi.quit()
snow::stopCluster(cl)

