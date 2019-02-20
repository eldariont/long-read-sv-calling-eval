library(tidyverse)

args = commandArgs(trailingOnly=TRUE)

res <- read_tsv(args[1], col_names = c("tool", "sample", "min", "metric", "value"))

res %>%
    mutate(fraction = str_match(sample, ".*\\.subsampled\\.(.*)")[2]) %>%
    spread(key = metric, value = value) %>%
    ggplot(aes(precision, recall, color=sample)) +
      geom_point(size=1.5) +
      geom_line() +
      theme_bw()

ggsave(args[2], width=10, height=5)
