library(tidyverse)

summary <- read.csv("small_sum.csv")
rownames(summary) <- summary[,1]
summary <- summary[,-1]

### create some tables to export to sqlite3, where we restrict to only significant probes
### i.e difference in beta expression of 33% and a p-adjusted of less than 0.01

probes_cl1 <- as.data.frame((summary$c1_avg)/(summary$norm_avg))
colnames(probes_cl1) <- ("diff")
rownames(probes_cl1) <- rownames(summary)
probes_cl1$pval <- summary$pval_norm_vs_c1

probes_cl2 <- as.data.frame((summary$c2_avg)/(summary$norm_avg))
colnames(probes_cl2) <- ("diff")
rownames(probes_cl2) <- rownames(summary)
probes_cl2$pval <- summary$pval_norm_vs_c2

probes_cl3 <- as.data.frame((summary$c3_avg)/(summary$norm_avg))
colnames(probes_cl3) <- ("diff")
rownames(probes_cl3) <- rownames(summary)
probes_cl3$pval <- summary$pval_norm_vs_c3

# write.csv(probes_cl1,"probe_cl1.csv")
# write.csv(probes_cl2,"probe_cl2.csv")
# write.csv(probes_cl3,"probe_cl3.csv")

###################   sig_probe.sh   ######################

cl1B <- as.data.frame(rownames(probes_cl1))
cl2B <- as.data.frame(rownames(probes_cl2))
cl3B <- as.data.frame(rownames(probes_cl3))

colnames(cl1B) <- c("probes")
colnames(cl2B) <- c("probes")
colnames(cl3B) <- c("probes")

sig_probes <- rbind(cl1B, cl2B, cl3B)
sig_probes <- as.data.frame(sig_probes[!duplicated(sig_probes), ])
colnames(sig_probes) <- c("probes")

###################   hyper hypo analysis  ###########3###########

hyper_cl1B <- probes_cl1[probes_cl1$diff > 1.03,]
hypo_cl1B <- probes_cl1[probes_cl1$diff < 0.97,]
hyper_cl2B <- probes_cl2[probes_cl2$diff > 1.03,]
hypo_cl2B <- probes_cl2[probes_cl2$diff < 0.97,]
hyper_cl3B <- probes_cl3[probes_cl3$diff > 1.03,]
hypo_cl3B <- probes_cl3[probes_cl3$diff < 0.97,]

hyper_cl1b <- as.data.frame(rownames(hyper_cl1B))
hypo_cl1b <- as.data.frame(rownames(hypo_cl1B))
hyper_cl2b <- as.data.frame(rownames(hyper_cl2B))
hypo_cl2b <- as.data.frame(rownames(hypo_cl2B))
hyper_cl3b <- as.data.frame(rownames(hyper_cl3B))
hypo_cl3b <- as.data.frame(rownames(hypo_cl3B))

colnames(hyper_cl1b) <- c("probes")
colnames(hypo_cl1b) <- c("probes")
colnames(hyper_cl2b) <- c("probes")
colnames(hypo_cl2b) <- c("probes")
colnames(hyper_cl3b) <- c("probes")
colnames(hypo_cl3b) <- c("probes")

###################   probe_to_gene.sh  ######################

manifest <- read_csv("headless.csv")

cl_temp <- merge(x = hyper_cl1b, y = manifest, by.x="probes", by.y="IlmnID", all = FALSE)
cl_temp <- as.data.frame(cl_temp$UCSC_RefGene_Name)
colnames(cl_temp) <- c("genes")
cl_genes <- cl_temp %>%
  separate_rows(genes, sep=";")
hyper_cl1_dm_genes <- unique(cl_genes)

cl_temp <- merge(x = hypo_cl1b, y = manifest, by.x="probes", by.y="IlmnID", all = FALSE)
cl_temp <- as.data.frame(cl_temp$UCSC_RefGene_Name)
colnames(cl_temp) <- c("genes")
cl_genes <- cl_temp %>%
  separate_rows(genes, sep=";")
hypo_cl1_dm_genes <- unique(cl_genes)

hyper_hypo_overlap_cl1 <- merge(x = hyper_cl1_dm_genes, y = hypo_cl1_dm_genes, by.x="genes", by.y="genes", all = FALSE)

cl_temp <- merge(x = hyper_cl2b, y = manifest, by.x="probes", by.y="IlmnID", all = FALSE)
cl_temp <- as.data.frame(cl_temp$UCSC_RefGene_Name)
colnames(cl_temp) <- c("genes")
cl_genes <- cl_temp %>%
  separate_rows(genes, sep=";")
hyper_cl2_dm_genes <- unique(cl_genes)

cl_temp <- merge(x = hypo_cl2b, y = manifest, by.x="probes", by.y="IlmnID", all = FALSE)
cl_temp <- as.data.frame(cl_temp$UCSC_RefGene_Name)
colnames(cl_temp) <- c("genes")
cl_genes <- cl_temp %>%
  separate_rows(genes, sep=";")
hypo_cl2_dm_genes <- unique(cl_genes)

hyper_hypo_overlap_cl2 <- merge(x = hyper_cl2_dm_genes, y = hypo_cl2_dm_genes, by.x="genes", by.y="genes", all = FALSE)

cl_temp <- merge(x = hyper_cl3b, y = manifest, by.x="probes", by.y="IlmnID", all = FALSE)
cl_temp <- as.data.frame(cl_temp$UCSC_RefGene_Name)
colnames(cl_temp) <- c("genes")
cl_genes <- cl_temp %>%
  separate_rows(genes, sep=";")
hyper_cl3_dm_genes <- unique(cl_genes)

cl_temp <- merge(x = hypo_cl3b, y = manifest, by.x="probes", by.y="IlmnID", all = FALSE)
cl_temp <- as.data.frame(cl_temp$UCSC_RefGene_Name)
colnames(cl_temp) <- c("genes")
cl_genes <- cl_temp %>%
  separate_rows(genes, sep=";")
hypo_cl3_dm_genes <- unique(cl_genes)

hyper_hypo_overlap_cl3 <- merge(x = hyper_cl3_dm_genes, y = hypo_cl3_dm_genes, by.x="genes", by.y="genes", all = FALSE)

######  make some files for GO enrichment analysis  ######

write.csv(hyper_cl1_dm_genes,"cl1_dm_genes_hyper.csv")
write.csv(hypo_cl1_dm_genes,"cl1_dm_genes_hypo.csv")
write.csv(hyper_cl2_dm_genes,"cl2_dm_genes_hyper.csv")
write.csv(hypo_cl2_dm_genes,"cl2_dm_genes_hypo.csv")
write.csv(hyper_cl3_dm_genes,"cl3_dm_genes_hyper.csv")
write.csv(hypo_cl3_dm_genes,"cl3_dm_genes_hypo.csv")

