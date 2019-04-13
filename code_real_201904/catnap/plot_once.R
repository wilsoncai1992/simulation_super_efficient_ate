load("./df_results.rda")
library(ggplot2)
gg <- ggplot(df_results, aes(x = as.factor(method), y = psi)) +
  geom_point() +
  geom_errorbar(aes(ymin=ci_lower, ymax=ci_upper, width=0.2)) +
  geom_hline(yintercept = 0, lty = 2) +
  facet_wrap(. ~ A, ncol = 4, scales = 'free') +
  theme_bw()
ggsave(gg, filename = 'once.png', width = 8, height = 8)