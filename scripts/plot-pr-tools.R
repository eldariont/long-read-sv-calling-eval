library(tidyverse)

args = commandArgs(trailingOnly=TRUE)

res <- read_tsv(args[1], col_names = c("tool", "sample", "min", "metric", "value"))

res %>%
    spread(key = metric, value = value) %>%
    filter(precision!=0 | recall!=0) %>%
    ggplot(aes(precision, recall, color=tool)) +
      geom_point(size=1.5) +
      geom_path() +
      coord_cartesian(xlim=c(0,1), ylim=c(0,1)) +
      theme_bw()

ggsave(args[2], width=6, height=5)
