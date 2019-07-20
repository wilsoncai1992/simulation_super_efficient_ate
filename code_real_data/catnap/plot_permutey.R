library(ggpubr)
# load("./result_permutey/permute_y.rda")
load("./result_npboot/npbootstrap.rda")
library(tidyverse)
# df_summary <- df_summary %>% filter(cnt >= 0.8 * max(cnt))
df_summary <- df_summary %>% mutate(method =
  recode(
    method,
    onestep_regular = "OS",
    onestep_reduced = "COS",
    tmle_regular = "TMLE",
    tmle_reduced = "CTMLE"
  )) %>% mutate(
    method = factor(method, levels = c("TMLE", "CTMLE", "OS", "COS"))
  )

gg_bias <- ggplot(
  df_summary %>% filter(method != "plugin"),
  aes(y = abs(bias), x = method)
) +
  geom_boxplot() +
  geom_hline(yintercept = 0, lty = 3) +
  ylab("Absolute bias (log-scale)") +
  scale_y_log10() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  rremove("x.text") +
  rremove("xlab")

gg_variance2 <- ggplot(
    df_summary %>% filter(method != "plugin"),
    aes(y = variance, x = method)
  ) +
  geom_boxplot() +
  ylab("Variance (log-scale)") +
  scale_y_log10(limits = c(NA, 1e1)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  rremove("xlab")
gg_variance <- gg_variance2 +
  rremove("x.text")

gg_mse <- ggplot(
    df_summary %>% filter(method != "plugin"),
    aes(y = mse, x = method)
  ) +
  ylab("MSE (log-scale)") +
  geom_boxplot() +
  theme_bw() +
  scale_y_log10(limits = c(NA, 1e1)) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  rremove("xlab")
gg_coverage <- ggplot(
    df_summary %>% filter(method != "plugin"),
    aes(y = coverage, x = method)
  ) +
  geom_boxplot() +
  ylab("Coverage") +
  theme_bw() +
  geom_hline(yintercept = .95, lty = 3) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  rremove("xlab")


gg1 <- ggarrange(gg_bias, gg_variance, gg_mse, gg_coverage, align = 'v', nrow = 2, ncol = 2, heights = c(1, 1.3))
ggsave(gg1, filename = 'permutey.pdf', width = 5, height = 4)

df1 <- df_summary %>% filter(method == "OS")
df2 <- df_summary %>% filter(method == "COS")
df3 <- df_summary %>% filter(method == "TMLE")
df4 <- df_summary %>% filter(method == "CTMLE")
# how many times onestep_reduced beat onestep_regular
# table(df1$mse >= df2$mse)
# table(df1$variance >= df2$variance)
# how many times tmle_reduced beat tmle_regular
# table(df3$mse >= df4$mse)
# table(df3$variance >= df4$variance)

ratio1 <- data.frame(var_ratio = df2$variance / df1$variance, method = 'OS', A = df2$A)
ratio2 <- data.frame(var_ratio = df4$variance / df3$variance, method = 'TMLE', A = df4$A)
ratios <- dplyr::bind_rows(ratio1, ratio2)

ggplot(ratios, aes(x = method, y = var_ratio)) +
  geom_boxplot() +
  scale_y_log10() +
  ylab("Relative efficiency (log-scale)") +
  theme_bw()

load("./smoothing/df_results_05.rda")
ratios <- left_join(ratios, df_positivity, by = "A")
gg_ratio <- ggplot(ratios, aes(lty = method, pch = method, y = var_ratio, x = positivity_score)) +
  geom_point() +
  geom_smooth(color = 'black') +
  scale_y_log10() +
  geom_hline(yintercept = 1, lty = 1) +
  ylab("Relative efficiency (log-scale)") +
  xlab(expression(P[n] (bar(G)[n] %in% (list(0.05, 0.95))))) +
  theme_bw() +
  theme(legend.position = "bottom") +
  guides(lty = guide_legend(title = "Method")) +
  guides(pch = guide_legend(title = "Method"))
ggsave(
  # gg_ratio,
  # filename = '0.05.pdf', width = 3, height = 3
  ggarrange(gg_variance2, gg_ratio, nrow = 1, ncol = 2, common.legend = TRUE, legend = "bottom", labels = "AUTO"),
  filename = '0.05.pdf', width = 7, height = 4
)
