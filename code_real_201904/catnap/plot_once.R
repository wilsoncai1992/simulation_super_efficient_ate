library(ggplot2)
gg <- ggplot(df_results, aes(x = as.factor(method), y = psi)) + 
  geom_point() + 
  geom_errorbar(aes(ymin=ci_lower, ymax=ci_upper, width=0.2)) +
  facet_wrap(. ~ A, ncol = 4, scales = 'free')
ggsave(gg, filename = 'once.png', width = 6, height = 6)