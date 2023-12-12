library(tidyverse)

summary <- read.csv("small_sum.csv")
rownames(summary) <- summary[,1]
summary <- summary[,-1]

### create some tables to export to sqlite3, where we restrict to only significant probes
### i.e difference in beta expression of 33% and a p-adjusted of less than 0.01

probes_cl1 <- data.frame(rownames(summary))
colnames(probes_cl1) <- c("probes")
probes_cl1$diff <- (summary$c1_avg)/(summary$norm_avg)
probes_cl1$pval <- summary$pval_norm_vs_c1

probes_cl2 <- data.frame(rownames(summary))
colnames(probes_cl2) <- c("probes")
probes_cl2$diff <- (summary$c2_avg)/(summary$norm_avg)
probes_cl2$pval <- summary$pval_norm_vs_c2

probes_cl3 <- data.frame(rownames(summary))
colnames(probes_cl3) <- c("probes")
probes_cl3$diff <- (summary$c3_avg)/(summary$norm_avg)
probes_cl3$pval <- summary$pval_norm_vs_c3

###################   probe_to_gene.sh  ######################

manifest <- read_csv("headless.csv")

cl_temp <- merge(x = probes_cl1, y = manifest, by.x="probes", by.y="IlmnID", all = FALSE)
cl_temp <- as.data.frame(cl_temp[, c('probes', 'UCSC_RefGene_Name', 'diff', 'pval')])
temp <- cl_temp %>%               
  separate_rows(UCSC_RefGene_Name, sep=";")
temp <- as.data.frame(unique(temp))
cl1_dm_genes <- temp[order(-temp$diff),]

cl_temp <- merge(x = probes_cl2, y = manifest, by.x="probes", by.y="IlmnID", all = FALSE)
cl_temp <- as.data.frame(cl_temp[, c('probes', 'UCSC_RefGene_Name', 'diff', 'pval')])
temp <- cl_temp %>%               
  separate_rows(UCSC_RefGene_Name, sep=";")
temp <- as.data.frame(unique(temp))
cl2_dm_genes <- temp[order(-temp$diff),]

cl_temp <- merge(x = probes_cl3, y = manifest, by.x="probes", by.y="IlmnID", all = FALSE)
cl_temp <- as.data.frame(cl_temp[, c('probes', 'UCSC_RefGene_Name', 'diff', 'pval')])
temp <- cl_temp %>%               
  separate_rows(UCSC_RefGene_Name, sep=";")
temp <- as.data.frame(unique(temp))
cl3_dm_genes <- temp[order(-temp$diff),]

write.csv(cl1_dm_genes, "cl1_genes.csv")
write.csv(cl2_dm_genes, "cl2_genes.csv")
write.csv(cl3_dm_genes, "cl3_genes.csv")

