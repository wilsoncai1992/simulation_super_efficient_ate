load("./output/df_mc_result.rda")
library(ggplot2)
library(ggpubr)
library(dplyr)

df_table <- df_mc_result %>% mutate(data_method = paste(dat_name, method, sep = "_")) %>% ungroup()
df_table <- df_table %>%
  filter(method %in% c(
    "oat_cv_variance",
    "oat_oracle_ci",
    "tmle_cv_variance",
    "tmle_oracle_ci",
    "onestep_oat_cv_variance",
    "onestep_oat_oracle_ci",
    "onestep_cv_variance",
    "onestep_oracle_ci"
  ))

df_bias <- df_table %>%
  select(dat_name, method, bias) %>%
  tidyr::spread(dat_name, bias) %>%
  mutate(metric = "bias")
df_variance <- df_table %>%
  select(dat_name, method, variance) %>%
  tidyr::spread(dat_name, variance) %>%
  mutate(metric = "variance")
df_mse <- df_table %>%
  select(dat_name, method, mse) %>%
  tidyr::spread(dat_name, mse) %>%
  mutate(metric = "mse")
df_coverage <- df_table %>%
  select(dat_name, method, coverage) %>%
  tidyr::spread(dat_name, coverage) %>%
  mutate(metric = "coverage")
# df_bias <- df_table %>%
#   select(data_method, bias) %>%
#   tidyr::spread(data_method, bias)
# df_bias <- as.data.frame(df_bias)
# rownames(df_bias) <- "bias"
# df_variance <- df_table %>%
#   select(data_method, variance) %>%
#   tidyr::spread(data_method, variance)
# df_variance <- as.data.frame(df_variance)
# rownames(df_variance) <- "variance"
# df_mse <- df_table %>%
#   select(data_method, mse) %>%
#   tidyr::spread(data_method, mse)
# df_mse <- as.data.frame(df_mse)
# rownames(df_mse) <- "mse"
# df_coverage <- df_table %>%
#   select(data_method, coverage) %>%
#   tidyr::spread(data_method, coverage)
# df_coverage <- as.data.frame(df_coverage)
# rownames(df_coverage) <- "coverage"

df_output <- rbind(df_bias, df_variance, df_mse, df_coverage)

# load("./output/df_no_permute.rda")
# df_once <- df_no_permute %>%
#   mutate(upper = Psi + 1.96 * sd, lower = Psi - 1.96 * sd)
# df_once <- df_once %>%
#   filter(method %in% c(
#     "oat_cv_variance",
#     "tmle_cv_variance",
#     "onestep_oat_cv_variance",
#     "onestep_cv_variance"
#   )) %>% select(dat_name, method, Psi, upper, lower)

save(
  df_output,
  # df_once,
  file = "./output/report.rda"
)
