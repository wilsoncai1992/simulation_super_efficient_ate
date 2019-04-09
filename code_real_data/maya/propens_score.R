library(dplyr)
load("./data/schwab.RData")
df_schwab <- data %>%
  mutate(
    Y = is.na(hiv_0),
    sex_0 = as.numeric(sex_0),
    formal_hi_occup_0 = as.numeric(formal_hi_occup_0),
    informal_hi_occup_0 = as.numeric(informal_hi_occup_0),
    informal_low_occup_0 = as.numeric(informal_low_occup_0),
    jobless_0 = as.numeric(jobless_0),
    edu_primary_0 = as.numeric(edu_primary_0),
    edu_secondary_plus_0 = as.numeric(edu_secondary_plus_0),
    mobile_0 = as.numeric(mobile_0),
    student_0 = as.numeric(student_0),
    # factor to numeric integer
    wealth_0 = as.numeric(as.character(wealth_0)),
    marital_0 = as.numeric(as.character(marital_0))
  ) %>%
  select(-hiv_0) %>%
  sample_frac(0.1) %>%
  filter(complete.cases(.))
library(hal9001)

# library(psych)
# pairs.panels(
#   df_schwab,
#   method = "pearson",
#   density = TRUE,
#   ellipses = TRUE
# )


W_names <- c(
  "sex_0",
  "age_0",
  "marital_0",
  # "formal_hi_occup_0",
  "informal_hi_occup_0",
  "informal_low_occup_0",
  "jobless_0",
  "edu_primary_0",
  "edu_secondary_plus_0",
  "wealth_0",
  "mobile_0",
  "student_0"
)
set_job <- c(
  # "formal_hi_occup_0",
  "informal_hi_occup_0",
  "informal_low_occup_0",
  "jobless_0"
)
set_edu <- c(
  "edu_primary_0",
  "edu_secondary_plus_0"
)

library(future)
library(listenv)
plan(multisession, workers = 8)
# plan(sequential)
df_propensity_scores <- listenv()

j <- 1
for (A_one in W_names) {
  df_propensity_scores[[j]] %<-% {
    W_remain <- W_names[-which(A_one == W_names)]
    if (A_one %in% set_job) {
      # if A is regarding job; remove any job correlated covariate from W
      W_remain <- setdiff(W_remain, set_job)
    }
    if (A_one %in% set_edu) {
      # if A is regarding edu; remove any edu correlated covariate from W
      W_remain <- setdiff(W_remain, set_edu)
    }
    df_W <- df_schwab[, W_remain]
    df_A <- df_schwab[, A_one]
    if (length(unique(df_A)) != 2) {
      # the outcome is not binary
      return(NULL)
    } else {
      g_fit <- hal9001::fit_hal(
        X = df_W,
        Y = as.numeric(df_A),
        # degrees = NULL,
        degrees = 2,
        fit_type = "glmnet",
        n_folds = 3,
        family = "binomial",
        return_lasso = TRUE,
        return_x_basis = TRUE,
        cv_select = TRUE,
        yolo = FALSE
      )
      g_hat <- predict(g_fit, new_data = df_W)
      data.frame(g_hat = g_hat, A_name = A_one)
    }
  }
  j <- j + 1
}
df_propensity_scores <- do.call(rbind, as.list(df_propensity_scores))
table(df_propensity_scores$A_name)
library(ggplot2)
gg <- ggplot(df_propensity_scores, aes(x = g_hat)) +
  geom_density() +
  facet_wrap(. ~ A_name, ncol = 3, scales = "free") +
  theme_bw()
ggsave(filename = "./g_fit.png", plot = gg, width = 6, height = 6)
