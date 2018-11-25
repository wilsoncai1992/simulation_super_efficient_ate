load("./output/df_mc_result.rda")
library(ggplot2)

ggplot(data = df_mc_result, aes(x = n, y = bias, color = method)) +
  geom_line() +
  scale_y_log10() +
  ylim(c(-1e-1, 1e-1))

ggplot(data = df_mc_result, aes(x = n, y = variance, color = method)) +
  geom_line() +
  scale_y_log10() +
  ylim(c(0, 1e-1))

ggplot(data = df_mc_result, aes(x = n, y = mse, color = method)) +
  geom_line() +
  scale_y_log10() +
  ylim(c(0, 1e-1))

ggplot(data = df_mc_result, aes(x = n, y = coverage, color = method)) +
  geom_line() +
  geom_hline(yintercept = .95, lty = 2)

library(dplyr)
yo <- df_simulation_result %>% filter(method == "tmle" & n == 1e3)
hist(yo$bias, 1e2)
hist(yo$mse)

yo <- df_simulation_result %>% filter(method == "onestep" & n == 1e3)
hist(log10(yo$bias), 1e2)
summary(yo$bias)
