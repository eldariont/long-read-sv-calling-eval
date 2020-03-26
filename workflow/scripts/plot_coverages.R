library(tidyverse)
library(scales)

args = commandArgs(trailingOnly=TRUE)

res <- read_tsv(args[1], col_names = c("caller", "mapper", "coverage", "zygosity", "type", "score", "metric", "value"))
res$type = factor(res$type, levels=c('del','inv','dup','ins'), labels=c('Deletions','Inversions','Tandem Dupl.','Insertions'))
res$caller = factor(res$caller, levels=c('pbsv', 'Sniffles', 'SVIM'), labels=c('pbsv', 'Sniffles', 'SVIM'))
res$zygosity = factor(res$zygosity, levels=c('homozygous', 'heterozygous'), labels=c('homozygous', 'heterozygous'))

res %>%
    filter(metric %in% c("recall", "precision")) %>%
    pivot_wider(names_from=metric, values_from=value) %>%
    filter(recall!=0 | precision!=0) %>%
    mutate(precision = 100*precision, recall = 100*recall) %>%
    mutate(f1 = 2*precision*recall/(precision+recall)) %>%
    group_by(caller, mapper, coverage, zygosity, type) %>%
    summarise(bestf1 = max(f1, na.rm = TRUE)) %>%
    ggplot(aes(coverage, bestf1, color=caller, pch=caller, linetype=zygosity)) +
      geom_path(, alpha=0.7) +
      geom_point(size=1.0, alpha=0.7) +
      scale_shape_manual(values=c(15,16,17)) +
      scale_color_manual(values=c("deepskyblue3", "goldenrod2", "firebrick2")) +
      facet_wrap(~type) +
      labs(y = "Precision", x = "Recall", color = "Tool", pch = "Tool", linetype = "Zygosity") +
      lims(x=c(0,60), y=c(0,100)) +
      theme_bw() +
      theme(panel.spacing = unit(0.75, "lines")) +
      theme(text = element_text(size=14), axis.text.x = element_text(size=9), axis.text.y = element_text(size=9))

ggsave(args[2], width=5, height=4)
