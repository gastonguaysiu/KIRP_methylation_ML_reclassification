library("tidyverse")

em_genes <- read.csv("kegg_david.csv")
colnames(em_genes) <- c("category", "term_description", "count", "percent", "p_val", "Benjamini")

ggplot(data=em_genes, aes(x = fct_reorder(term_description, percent), y=percent)) + 
  geom_col(aes(fill=p_val))  +
  ylim(0, 15) +
  labs(title= "Enrichment",
       caption = "percent = number of associated genes on list found in David / all genes listed found in David") + 
  theme(
    axis.title.y = element_text(size = 8),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    legend.title = element_text(size = 10), 
    legend.text  = element_text(size = 8),
    legend.key.size = unit(0.5, "lines")) +
  coord_flip()
