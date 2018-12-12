load("./output/df_mc_result.rda")
library(ggplot2)
library(ggpubr)
library(dplyr)

df_plot_oat1 <- df_mc_result %>%
  filter(method %in% c("oat_oracle_ci", "tmle"))
df_plot_oat1$method <- as.character(df_plot_oat1$method)
df_plot_oat1$method[df_plot_oat1$method == "oat_oracle_ci"] <- "ctmle"
gg1 <- ggplot(
  data = df_plot_oat1, aes(x = n, y = bias, shape = method, lty = method)
) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0, lty = 3) +
  facet_grid(. ~ iv_beta) +
  theme_bw()
gg2 <- ggplot(
  data = df_plot_oat1, aes(x = n, y = variance, shape = method, lty = method)
) +
  geom_line() +
  geom_point() +
  facet_grid(. ~ iv_beta) +
  theme_bw()

df_tmle <- df_plot_oat1 %>% filter(method %in% c("tmle"))
df_plot <- dplyr::left_join(df_plot_oat1, df_tmle, c("n", "iv_beta"))
df_plot$re <- df_plot$mse.x / df_plot$mse.y
df_plot$method <- df_plot$method.x
gg3 <- ggplot(
  data = df_plot, aes(x = n, y = re, shape = method, lty = method)
) +
  geom_line() +
  geom_point() +
  facet_grid(. ~ iv_beta) +
  theme_bw()
gg_panel1 <- ggarrange(
  gg1,
  gg2,
  gg3,
  nrow = 3,
  labels = "AUTO",
  common.legend = TRUE,
  legend = "none"
)

# df_sdOracle <- df_simulation_result %>%
#   group_by(n, iv_beta, method) %>%
#   summarise(sd_oracle = sd(bias))
# df_bias_std <- dplyr::left_join(
#   df_simulation_result, df_sdOracle, c("n", "iv_beta", "method")
# )
# df_bias_std$bias_std <- df_bias_std$bias / df_bias_std$sd_oracle
df_plot_density <- df_simulation_result %>%
  filter(method %in% c("oat_cv_variance", "tmle"))
df_plot_density$method <- as.character(df_plot_density$method)
df_plot_density$method[df_plot_density$method == "oat_cv_variance"] <- "ctmle"

gg4 <- ggplot(
  data = df_plot_density, aes(x = bias * sqrt(n), lty = method)
) +
  geom_density(alpha = 1) +
  geom_vline(xintercept = 0, lty = 3) +
  xlim(c(-1, 1) * 10) +
  facet_grid(n ~ iv_beta) +
  theme_bw() + theme(legend.position = "none")
gg_panel2 <- ggarrange(
  gg_panel1,
  gg4,
  ncol = 2,
  labels = c("", "D"),
  common.legend = TRUE,
  legend = "bottom"
)
ggsave(gg_panel2, filename = "./output/tmle_panel.pdf", width = 12, height = 6.75)

gg <- ggplot(
  data = df_plot_oat1, aes(x = n, y = coverage, shape = method, lty = method)
) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = .95, lty = 2) +
  facet_grid(. ~ iv_beta) +
  theme_bw() +
  theme(legend.position = "bottom")
ggsave(gg, filename = "./output/tmle_coverage.pdf", width = 4, height = 4)

# =============================================================================
df_plot_oat1 <- df_mc_result %>%
  filter(method %in% c("onestep_oat_oracle_ci", "onestep"))
df_plot_oat1$method <- as.character(df_plot_oat1$method)
df_plot_oat1$method[df_plot_oat1$method == "onestep_oat_oracle_ci"] <- "c-onestep"


gg1 <- ggplot(
  data = df_plot_oat1, aes(x = n, y = bias, shape = method, lty = method)
) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0, lty = 3) +
  facet_grid(. ~ iv_beta) +
  theme_bw()
gg2 <- ggplot(
  data = df_plot_oat1, aes(x = n, y = variance, shape = method, lty = method)
) +
  geom_line() +
  geom_point() +
  facet_grid(. ~ iv_beta) +
  theme_bw()

df_tmle <- df_plot_oat1 %>% filter(method %in% c("onestep"))
df_plot <- dplyr::left_join(df_plot_oat1, df_tmle, c("n", "iv_beta"))
df_plot$re <- df_plot$mse.x / df_plot$mse.y
df_plot$method <- df_plot$method.x
gg3 <- ggplot(
  data = df_plot, aes(x = n, y = re, shape = method, lty = method)
) +
  geom_line() +
  geom_point() +
  facet_grid(. ~ iv_beta) +
  theme_bw()
gg_panel1 <- ggarrange(
  gg1,
  gg2,
  gg3,
  nrow = 3,
  labels = "AUTO",
  common.legend = TRUE,
  legend = "none"
)

df_plot_density <- df_simulation_result %>%
  filter(method %in% c("onestep_oat_oracle_ci", "onestep"))
df_plot_density$method <- as.character(df_plot_density$method)
df_plot_density$method[df_plot_density$method == "onestep_oat_oracle_ci"] <- "c-onestep"

gg4 <- ggplot(
  data = df_plot_density, aes(x = bias * sqrt(n), lty = method)
) +
  geom_density(alpha = 1) +
  geom_vline(xintercept = 0, lty = 3) +
  xlim(c(-1, 1) * 10) +
  facet_grid(n ~ iv_beta) +
  theme_bw() + theme(legend.position = "none")
gg_panel2 <- ggarrange(
  gg_panel1,
  gg4,
  ncol = 2,
  labels = c("", "D"),
  common.legend = TRUE,
  legend = "bottom"
)
ggsave(gg_panel2, filename = "./output/onestep_panel.pdf", width = 12, height = 6.75)

gg <- ggplot(
  data = df_plot_oat1, aes(x = n, y = coverage, shape = method, lty = method)
) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = .95, lty = 2) +
  facet_grid(. ~ iv_beta) +
  theme_bw() +
  theme(legend.position = "bottom")
ggsave(gg, filename = "./output/onestep_coverage.pdf", width = 4, height = 4)
