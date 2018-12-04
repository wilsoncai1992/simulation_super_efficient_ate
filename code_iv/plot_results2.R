load("./output/df_mc_result.rda")
library(ggplot2)
library(dplyr)

gg <- ggplot(data = df_mc_result, aes(x = n, y = bias, color = method)) +
  geom_line() +
  geom_point() +
  facet_grid(. ~ iv_beta)
  # ylim(c(-1e-1, 1e-1))
ggsave(filename = "./output/bias.png", plot = gg, width = 6, height = 6)

gg <- ggplot(data = df_mc_result, aes(x = n, y = variance, color = method)) +
  geom_line() +
  geom_point() +
  facet_grid(. ~ iv_beta)
  # ylim(c(0, 1e-1))
ggsave(filename = "./output/variance.png", plot = gg, width = 6, height = 6)

gg <- ggplot(data = df_mc_result, aes(x = n, y = mse, color = method)) +
  geom_line() +
  geom_point() +
  facet_grid(. ~ iv_beta)
  # ylim(c(0, 1e-1))
ggsave(filename = "./output/mse.png", plot = gg, width = 6, height = 6)

gg <- ggplot(
  data = df_mc_result %>% filter(!grepl(pattern = 'plugin', method)),
  aes(x = n, y = coverage, color = method)
) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = .95, lty = 2) +
  facet_grid(. ~ iv_beta)
ggsave(filename = "./output/coverage.png", plot = gg, width = 6, height = 6)

gg <- ggplot(
  data = df_mc_result %>% filter(!grepl(pattern = 'plugin', method)),
  aes(x = as.factor(n), y = sd, color = method)
) +
  geom_boxplot() +
  facet_grid(. ~ iv_beta)
ggsave(filename = "./output/sd.png", plot = gg, width = 6, height = 6)
