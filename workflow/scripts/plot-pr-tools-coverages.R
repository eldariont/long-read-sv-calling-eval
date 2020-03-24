library(tidyverse)

args = commandArgs(trailingOnly=TRUE)

res <- read_tsv(args[1], col_names = c("tool", "sample", "min", "metric", "value"))
cov <- read_tsv(args[2], col_names = c("sample", "coverage"))
cov$coverage <- as.factor(cov$coverage)

res %>%
    inner_join(cov, by = "sample") %>%
    spread(key = metric, value = value) %>%
    filter(precision != 0 | recall != 0) %>%
    ggplot(aes(precision, recall, color=coverage)) +
      geom_point(size=1.5) +
      geom_path(aes(linetype=coverage)) +
      coord_cartesian(xlim=c(0,1), ylim=c(0,1)) +
      theme_bw()

ggsave(args[3], width=6, height=5)
