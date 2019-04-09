library(dplyr)
library(glmnet)
library(future)
library(listenv)
load("./all_data.RData")
A_names <- names(all_data)

plan(multisession, workers = 8)
# plan(sequential)
df_propensity_scores <- listenv()

j <- 1
for (A_name in A_names) {
  df_propensity_scores[[j]] %<-% {
    data_train <- all_data[[A_name]]
    # all W are integer type
    W <- as.matrix(data_train$W)
    A <- data_train$A
    is_not_missing <- complete.cases(W)
    W <- W[is_not_missing, ]
    A <- A[is_not_missing]
    lasso_fit <- cv.glmnet(
      x = W, y = A, nfolds = 5, family = "binomial", standardize = FALSE
    )
    # plot(lasso_fit)
    g_hat <- predict(lasso_fit, newx = W, s = "lambda.min", type = "response")
    # g_hat <- predict(lasso_fit, newx = W, s = "lambda.1se", type = "response")
    beta_hat <- lasso_fit$glmnet.fit$beta[
      , lasso_fit$glmnet.fit$lambda == lasso_fit$lambda.min
    ]
    # beta_hat <- lasso_fit$glmnet.fit$beta[
      # , lasso_fit$glmnet.fit$lambda == lasso_fit$lambda.1se
    # ]
    nonzero <- sum(beta_hat != 0)
    data.frame(g_hat = g_hat, A = A_name, nonzero = nonzero)
  }
  j <- j + 1
}
df_propensity_scores <- do.call(rbind, as.list(df_propensity_scores))
table(df_propensity_scores$A)
library(ggplot2)
# gg <- ggplot(df_propensity_scores, aes(x = g_hat)) +
gg <- ggplot(df_propensity_scores, aes(x = X1)) +
  geom_histogram() +
  facet_wrap(. ~ A, ncol = 8) +
  theme_bw()
df_beta_nonzero <- df_propensity_scores %>% 
  group_by(A) %>% 
  summarise(nonzero = max(nonzero))
ggsave(filename = "./g_fit.png", plot = gg, width = 10, height = 6)

write.csv(df_beta_nonzero, file = "df_beta_nonzero.csv")
