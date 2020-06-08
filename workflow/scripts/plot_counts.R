library(tidyverse)
library(scales)

args = commandArgs(trailingOnly=TRUE)

res <- read_tsv(args[1], col_names = c("size", "type", "count", "cumlength", "tool"))
res[res$type == "DUP", "type"] <- "DUP:TANDEM"
res[res$type == "DEL/INV", "type"] <- "Complex"
res[res$type == "DUP/INS", "type"] <- "Complex"
res[res$type == "INVDUP", "type"] <- "Complex"
res[res$type == "INV/INVDUP", "type"] <- "Complex"
res$size = factor(res$size, levels=c('tiny', 'small', 'medium', 'large', 'huge', 'all'), labels=c('tiny', 'small', 'medium', 'large', 'huge', 'all'))
res$type = factor(res$type, levels=c('DEL', 'INS', 'INV', 'DUP:TANDEM', 'DUP:INT', 'Complex'), labels=c('DEL', 'INS', 'INV', 'DUP:TANDEM', 'DUP:INT', 'Complex'))
res$tool = factor(res$tool, levels=c('pbsv', 'Sniffles', 'SVIM'), labels=c('pbsv', 'Sniffles', 'SVIM'))

res%>%
    ggplot(aes(tool, count, fill=type)) +
      geom_col() +
      scale_fill_brewer(palette="RdYlBu") +
      labs(y = "Count", x = "Tool", color = "SV class") +
      facet_grid(~size) +
      theme_bw() +
      theme(text = element_text(size=6), axis.text.x = element_text(size=6), axis.text.y = element_text(size=6))
ggsave(args[2], width=5, height=3)
res %>% write_tsv(args[3])

