load("./output/df_mc_result.rda")
library(ggplot2)

gg <- ggplot(data = df_mc_result, aes(x = n, y = bias, color = method)) +
  geom_line() +
  geom_point()
  # ylim(c(-1e-1, 1e-1))
ggsave(filename = "./output/bias.png", plot = gg, width = 6, height = 6)

gg <- ggplot(data = df_mc_result, aes(x = n, y = variance, color = method)) +
  geom_line() +
  geom_point()
  # ylim(c(0, 1e-1))
ggsave(filename = "./output/variance.png", plot = gg, width = 6, height = 6)

gg <- ggplot(data = df_mc_result, aes(x = n, y = mse, color = method)) +
  geom_line() +
  geom_point()
  # ylim(c(0, 1e-1))
ggsave(filename = "./output/mse.png", plot = gg, width = 6, height = 6)

gg <- ggplot(data = df_mc_result, aes(x = n, y = coverage, color = method)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = .95, lty = 2)
ggsave(filename = "./output/coverage.png", plot = gg, width = 6, height = 6)


gg <- ggplot(data = df_simulation_result, aes(x = as.factor(n), y = sd, color = method)) +
  geom_boxplot()
ggsave(filename = "./output/sd.png", plot = gg, width = 6, height = 6)
# library(dplyr)
# yo <- df_simulation_result %>% filter(method == "tmle" & n == 1e3)
# hist(yo$bias, 1e2)
# hist(yo$mse)
# 
# yo <- df_simulation_result %>% filter(method == "onestep" & n == 1e3)
# hist(log10(yo$bias), 1e2)
# summary(yo$bias)
