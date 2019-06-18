load("./all_data.RData")
source("./fit_one_treatment.R")
A_names <- names(all_data)

library(future)
library(listenv)
plan(multisession, workers = 8)
# plan(sequential)
df_results <- listenv()
j <- 1
for (A_name in A_names) {
  df_results[[j]] %<-% {
    data_train <- all_data[[A_name]]
    # all W are integer type
    W <- as.matrix(data_train$W)
    A <- data_train$A
    is_not_missing <- complete.cases(W)
    W <- W[is_not_missing, ]
    A <- A[is_not_missing]
    Y <- data_train$Y[is_not_missing]
    # fit_one_A(W = W, A = A, Y = Y, delta = 0.1)
    # fit_one_A(W = W, A = A, Y = Y, delta = 0.05)
    # fit_one_A(W = W, A = A, Y = Y, delta = 0.025)
    fit_one_A(W = W, A = A, Y = Y, delta = 0.01)
  }
  j <- j + 1
}

df_results <- do.call(rbind, as.list(df_results))
table(df_results$A)


library(tidyverse)
df_positivity <- df_results %>% group_by(A) %>% summarise(positivity_score = mean(positivity_score))
df_positivity$positivity_level <- cut(
  df_positivity$positivity_score,
  # breaks = c(quantile(df_positivity$positivity_score, probs = seq(0, 1, by = 0.25))),
  breaks = c(0, 7e-2, 1e-1, 1),
  include.lowest = TRUE
)
# save(df_results, df_positivity, file = "df_results.rda")
save(df_results, df_positivity, file = "df_results_01.rda")
