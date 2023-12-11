library(tidyverse)
library(stats)

island <- read.csv("CpG_Island.csv")
rownames(island) <- island[,1]
island <- island[,-1]

CpG_island_p <- fisher.test(island[,c(1,2)], simulate.p.value = TRUE, B = 10000)$p.value

gene <- read.csv("RefGene_Group.csv")
rownames(gene) <- gene[,1]
gene <- gene[,-1]

CpG_gene_p <- fisher.test(gene[,c(1,2)], simulate.p.value = TRUE, B = 10000)$p.value
