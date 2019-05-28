load("./all_data.RData")
source("./fit_one_treatment.R")
A_names <- names(all_data)

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
    fit_one_A(W = W, A = A, Y = Y)
  }
  j <- j + 1
}

df_results <- do.call(rbind, as.list(df_results))
table(df_results$A)

save(df_results, file = "df_results.rda")
