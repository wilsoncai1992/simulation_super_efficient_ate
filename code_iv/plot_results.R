load("./output/df_mc_result.rda")
library(ggplot2)
library(dplyr)

iv_beta_grid <- unique(df_mc_result$iv_beta)

do_one_iv_beta <- function(iv_beta){
  df_mc_result <- df_mc_result %>% filter(iv_beta == iv_beta)
  gg <- ggplot(data = df_mc_result, aes(x = n, y = bias, color = method)) +
    geom_line() +
    geom_point()
    # ylim(c(-1e-1, 1e-1))
  ggsave(filename = paste("./output/bias_", iv_beta, ".png"), plot = gg, width = 6, height = 6)

  gg <- ggplot(data = df_mc_result, aes(x = n, y = variance, color = method)) +
    geom_line() +
    geom_point()
    # ylim(c(0, 1e-1))
  ggsave(filename = paste("./output/variance_", iv_beta, ".png"), plot = gg, width = 6, height = 6)

  gg <- ggplot(data = df_mc_result, aes(x = n, y = mse, color = method)) +
    geom_line() +
    geom_point()
    # ylim(c(0, 1e-1))
  ggsave(filename = paste("./output/mse_", iv_beta, ".png"), plot = gg, width = 6, height = 6)

  gg <- ggplot(
    data = df_mc_result %>% filter(!grepl(pattern = 'plugin', method)),
    aes(x = n, y = coverage, color = method)
  ) +
    geom_line() +
    geom_point() +
    geom_hline(yintercept = .95, lty = 2)
  ggsave(filename = paste("./output/coverage_", iv_beta, ".png"), plot = gg, width = 6, height = 6)

  gg <- ggplot(
    data = df_mc_result %>% filter(!grepl(pattern = 'plugin', method)),
    aes(x = as.factor(n), y = sd, color = method)
  ) +
    geom_boxplot()
  ggsave(filename = paste("./output/sd_", iv_beta, ".png"), plot = gg, width = 6, height = 6)

}
