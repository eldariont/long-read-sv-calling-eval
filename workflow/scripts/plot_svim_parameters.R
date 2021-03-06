library(tidyverse)
library(scales)

args = commandArgs(trailingOnly=TRUE)

res <- read_tsv(args[1], col_names = c("parameters", "caller", "mapper", "subsample", "vcf", "score", "metric", "value"))

res %>%
    filter(metric %in% c("recall", "precision")) %>%
    pivot_wider(names_from=metric, values_from=value) %>%
    filter(recall!=0 | precision!=0) %>%
    mutate(precision = 100*precision, recall = 100*recall) %>%
    ggplot(aes(recall, precision, color=parameters, pch=parameters)) +
      geom_point(size=0.5) +
      geom_path() +
      facet_grid(subsample~vcf) +
      labs(y = "Precision", x = "Recall", color = "Tool", pch = "Tool") +
      lims(x=c(0,100), y=c(0,100)) +
      theme_bw() +
      theme(panel.spacing = unit(0.75, "lines")) +
      theme(text = element_text(size=14), axis.text.x = element_text(size=9), axis.text.y = element_text(size=9))

ggsave(args[2], width=20, height=12)
