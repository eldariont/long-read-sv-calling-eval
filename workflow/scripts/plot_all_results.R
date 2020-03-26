library(tidyverse)
library(scales)

args = commandArgs(trailingOnly=TRUE)

res <- read_tsv(args[1], col_names = c("caller", "mapper", "coverage", "zygosity", "type", "score", "metric", "value"))
res$type = factor(res$type, levels=c('del','inv','dup','ins'), labels=c('Deletions','Inversions','Tandem Dupl.','Insertions'))
res$caller = factor(res$caller, levels=c('Sniffles', 'SVIM', 'pbsv'), labels=c('Sniffles', 'SVIM', 'pbsv'))

res %>%
    filter(metric %in% c("recall", "precision")) %>%
    pivot_wider(names_from=metric, values_from=value) %>%
    filter(recall!=0 | precision!=0) %>%
    mutate(precision = 100*precision, recall = 100*recall) %>%
    ggplot(aes(recall, precision, color=caller, pch=caller)) +
      geom_point(size=0.5) +
      scale_shape_manual(values=c(15,16,17)) +
      scale_color_manual(values=c("deepskyblue", "dodgerblue4", "goldenrod1")) +
      geom_path() +
      facet_grid(coverage+zygosity~type) +
      labs(y = "Precision", x = "Recall", color = "Tool", pch = "Tool") +
      lims(x=c(0,100), y=c(0,100)) +
      theme_bw() +
      theme(panel.spacing = unit(0.75, "lines")) +
      theme(text = element_text(size=14), axis.text.x = element_text(size=9), axis.text.y = element_text(size=9))

ggsave(args[2], width=10, height=15)