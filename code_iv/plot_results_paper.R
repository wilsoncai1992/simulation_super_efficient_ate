load("./output/df_mc_result.rda")
library(ggplot2)
library(ggpubr)
library(dplyr)

df_mc_result <- df_mc_result %>% rename(Method = method)
df_simulation_result <- df_simulation_result %>% rename(Method = method)

get_df_plot <- function(df, method_to_plot, m1, m2) {
  df_out <- df %>% filter(Method %in% method_to_plot)
  df_out$Method <- as.character(df_out$Method)
  df_out$Method[df_out$Method == method_to_plot[1]] <- m1
  df_out$Method[df_out$Method == method_to_plot[2]] <- m2
  return(df_out)
}
df_plot_oat1 <- get_df_plot(
  df_mc_result, c("oat_oracle_ci", "tmle_oracle_ci"), "CTMLE", "TMLE"
)
df_plot_oat2 <- get_df_plot(
  df_mc_result, c("oat_cv_variance", "tmle_cv_variance"), "CTMLE", "TMLE"
)

gg1 <- ggplot(
  data = df_plot_oat1, aes(x = n, y = bias, shape = Method, lty = Method)
) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0, lty = 3) +
  ylab("Bias") +
  facet_grid(. ~ iv_beta) +
  theme_bw()
gg2 <- ggplot(
  data = df_plot_oat1, aes(x = n, y = variance, shape = Method, lty = Method)
) +
  geom_line() +
  geom_point() +
  ylab("Variance") +
  facet_grid(. ~ iv_beta) +
  theme_bw()

df_tmle <- df_plot_oat1 %>% filter(Method %in% c("TMLE"))
df_plot <- dplyr::left_join(df_plot_oat1, df_tmle, c("n", "iv_beta"))
df_plot$re <- df_plot$mse.x / df_plot$mse.y
df_plot$Method <- df_plot$Method.x
gg3 <- ggplot(
  data = df_plot, aes(x = n, y = re, shape = Method, lty = Method)
) +
  geom_line() +
  geom_point() +
  ylab("Relative efficiency") +
  facet_grid(. ~ iv_beta) +
  ylim(c(.1, 1.1)) +
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
#   group_by(n, iv_beta, Method) %>%
#   summarise(sd_oracle = sd(bias))
# df_bias_std <- dplyr::left_join(
#   df_simulation_result, df_sdOracle, c("n", "iv_beta", "Method")
# )
# df_bias_std$bias_std <- df_bias_std$bias / df_bias_std$sd_oracle
df_plot_density <- get_df_plot(
  df_simulation_result, c("oat_oracle_ci", "tmle_oracle_ci"), "CTMLE", "TMLE"
)
gg4 <- ggplot(
  data = df_plot_density, aes(x = bias * sqrt(n), lty = Method)
) +
  geom_density(alpha = 1) +
  geom_vline(xintercept = 0, lty = 3) +
  xlim(c(-1, 1) * 10) +
  xlab(expression(sqrt(n) * (Psi[n] - Psi[0]))) +
  ylab("Density") +
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
# ggsave(gg_panel2, filename = "./output/tmle_panel.png", width = 12, height = 6.75)

plot_coverage <- function(df) {
  ggplot(
    data = df, aes(x = n, y = coverage, shape = Method, lty = Method)
  ) +
    geom_line() +
    geom_point() +
    geom_hline(yintercept = .95, lty = 3) +
    ylab("Coverage") +
    ylim(c(.6, 1)) +
    facet_grid(. ~ iv_beta) +
    theme_bw()
}
gg01 <- plot_coverage(df_plot_oat1)
gg02 <- plot_coverage(df_plot_oat2)
gg_panel2 <- ggarrange(
  gg01,
  gg02,
  nrow = 2,
  labels = "AUTO",
  common.legend = TRUE,
  legend = "bottom"
)
ggsave(gg_panel2, filename = "./output/tmle_coverage.pdf", width = 4, height = 4)
# ggsave(gg_panel2, filename = "./output/tmle_coverage.png", width = 4, height = 4)

# =============================================================================
df_plot_oat1 <- get_df_plot(
  df_mc_result, c("onestep_oat_oracle_ci", "onestep_oracle_ci"), "C-onestep", "onestep"
)
df_plot_oat2 <- get_df_plot(
  df_mc_result, c("onestep_oat_cv_variance", "onestep_cv_variance"), "C-onestep", "onestep"
)

gg1 <- ggplot(
  data = df_plot_oat1, aes(x = n, y = bias, shape = Method, lty = Method)
) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0, lty = 3) +
  ylab("Bias") +
  facet_grid(. ~ iv_beta) +
  theme_bw()
gg2 <- ggplot(
  data = df_plot_oat1, aes(x = n, y = variance, shape = Method, lty = Method)
) +
  geom_line() +
  geom_point() +
  ylab("Variance") +
  facet_grid(. ~ iv_beta) +
  theme_bw()

df_tmle <- df_plot_oat1 %>% filter(Method %in% c("onestep"))
df_plot <- dplyr::left_join(df_plot_oat1, df_tmle, c("n", "iv_beta"))
df_plot$re <- df_plot$mse.x / df_plot$mse.y
df_plot$Method <- df_plot$Method.x
gg3 <- ggplot(
  data = df_plot, aes(x = n, y = re, shape = Method, lty = Method)
) +
  geom_line() +
  geom_point() +
  ylab("Relative efficiency") +
  facet_grid(. ~ iv_beta) +
  ylim(c(.1, 1.1)) +
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

df_plot_density <- get_df_plot(
  df_simulation_result, c("onestep_oat_oracle_ci", "onestep_oracle_ci"), "C-onestep", "onestep"
)
gg4 <- ggplot(
  data = df_plot_density, aes(x = bias * sqrt(n), lty = Method)
) +
  geom_density(alpha = 1) +
  geom_vline(xintercept = 0, lty = 3) +
  xlim(c(-1, 1) * 10) +
  xlab(expression(sqrt(n) * (Psi[n] - Psi[0]))) +
  ylab("Density") +
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
# ggsave(gg_panel2, filename = "./output/onestep_panel.png", width = 12, height = 6.75)

gg01 <- plot_coverage(df_plot_oat1)
gg02 <- plot_coverage(df_plot_oat2)
gg_panel2 <- ggarrange(
  gg01,
  gg02,
  nrow = 2,
  labels = "AUTO",
  common.legend = TRUE,
  legend = "bottom"
)
ggsave(gg_panel2, filename = "./output/onestep_coverage.pdf", width = 4, height = 4)
# ggsave(gg_panel2, filename = "./output/onestep_coverage.png", width = 4, height = 4)
