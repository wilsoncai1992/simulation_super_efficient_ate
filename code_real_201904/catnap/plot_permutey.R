# load("./permute_y.rda")
load("./npbootstrap.rda")
library(tidyverse)
gg_bias <- ggplot(
  df_summary %>% filter(method != "plugin"),
  aes(y = abs(bias), x = method)
) +
  geom_boxplot() +
  geom_hline(yintercept = 0, lty = 3) +
  # scale_y_sqrt() +
  scale_y_log10() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
gg_variance <- ggplot(
    df_summary %>% filter(method != "plugin"),
    aes(y = variance, x = method)
  ) +
  geom_boxplot() +
  scale_y_log10(limits = c(NA, 1e1)) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
gg_mse <- ggplot(
    df_summary %>% filter(method != "plugin"),
    aes(y = mse, x = method)
  ) +
  geom_boxplot() +
  scale_y_log10(limits = c(NA, 1e1)) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
gg_coverage <- ggplot(
    df_summary %>% filter(method != "plugin"),
    aes(y = coverage, x = method)
  ) +
  geom_boxplot() +
  geom_hline(yintercept = .95, lty = 3) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))


library(ggpubr)
gg1 <- ggarrange(gg_bias, gg_variance, gg_mse, align = 'v', nrow = 3)
gg2 <- ggarrange(gg1, gg_coverage, ncol = 2)
ggsave(gg2, filename = 'permutey.png', width = 4, height = 6)

df1 <- df_summary %>% filter(method == "onestep_regular")
df2 <- df_summary %>% filter(method == "onestep_reduced")
table(df1$mse <= df2$mse)

df1[df1$mse >= df2$mse, ]
