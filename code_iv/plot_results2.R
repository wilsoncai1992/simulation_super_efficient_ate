load("./output/df_mc_result.rda")
library(ggplot2)
library(dplyr)

gg <- ggplot(
  data = df_mc_result %>% filter(method %in% c("oat_oracle_ci", "oat_cv_variance", "tmle")),
  # data = df_mc_result %>% filter(!grepl(pattern = "plugin", method)),
  aes(x = n, y = bias, color = method)
) +
  geom_line() +
  geom_point() +
  # ylim(c(-.5, .5)) +
  facet_grid(. ~ iv_beta)
ggsave(filename = "./output/bias.png", plot = gg, width = 6, height = 6)

gg <- ggplot(
  data = df_mc_result %>% filter(method %in% c("oat_oracle_ci", "oat_cv_variance", "tmle")),
  # data = df_mc_result %>% filter(!grepl(pattern = "plugin", method)),
  aes(x = n, y = variance, color = method)
) +
  geom_line() +
  geom_point() +
  # ylim(c(0, 1)) +
  facet_grid(. ~ iv_beta)
ggsave(filename = "./output/variance.png", plot = gg, width = 6, height = 6)

df_oat <- df_mc_result %>% filter(method %in% c("oat_oracle_ci", "oat_cv_variance", "tmle"))
# df_oat <- df_mc_result %>% filter(!grepl(pattern = "plugin", method))
df_tmle <- df_mc_result %>% filter(method %in% c("tmle"))

df_plot <- dplyr::left_join(df_oat, df_tmle, c("n", "iv_beta"))
df_plot$re <- df_plot$mse.x / df_plot$mse.y
df_plot$method <- df_plot$method.x
gg <- ggplot(data = df_plot, aes(x = n, y = re, color = method)) +
  geom_line() +
  geom_point() +
  # ylim(c(0, 1)) +
  facet_grid(. ~ iv_beta)
ggsave(filename = "./output/re.png", plot = gg, width = 6, height = 6)

gg <- ggplot(
  data = df_mc_result %>% filter(method %in% c("oat_oracle_ci", "oat_cv_variance", "tmle")),
  # data = df_mc_result %>% filter(!grepl(pattern = "plugin", method)),
  aes(x = n, y = coverage, color = method)
) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = .95, lty = 2) +
  facet_grid(. ~ iv_beta)
ggsave(filename = "./output/coverage.png", plot = gg, width = 6, height = 6)


df_sdOracle <- df_simulation_result %>%
  group_by(n, iv_beta, method) %>%
  summarise(sd_oracle = sd(bias))
df_bias_std <- dplyr::left_join(
  df_simulation_result, df_sdOracle, c("n", "iv_beta", "method")
)
df_bias_std$bias_std <- df_bias_std$bias / df_bias_std$sd_oracle
gg <- ggplot(
  data = df_bias_std %>% filter(method %in% c("oat_oracle_ci", "oat_cv_variance", "tmle")),
  # data = df_bias_std %>% filter(!grepl(pattern = "plugin", method)),
  aes(x = bias_std, color = method)
) +
  stat_function(fun = function(x) dnorm(x), color = "black") +
  geom_density() +
  geom_vline(xintercept = 0, lty = 2) +
  facet_grid(n ~ iv_beta)
ggsave(filename = "./output/hist_bias_std.png", plot = gg, width = 6, height = 6)

gg <- ggplot(
  data = df_mc_result %>% filter(!grepl(pattern = "plugin", method)),
  data = df_mc_result %>% filter(method %in% c("oat_oracle_ci", "oat_cv_variance", "tmle")),
  # aes(x = as.factor(n), y = sd, color = method)
) +
  geom_boxplot() +
  ylim(c(0, 1)) +
  facet_grid(. ~ iv_beta)
ggsave(filename = "./output/sd.png", plot = gg, width = 6, height = 6)

write.csv(df_mc_result, "./output/table.csv")