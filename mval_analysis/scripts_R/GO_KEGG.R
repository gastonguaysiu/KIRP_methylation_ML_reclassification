
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)


gene_list <- read.csv("cl2_genes.csv")
genes_to_t <- rownames(gene_list)

GO_results <- enrichGO(gene = genes_to_t, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID", ont = "ALL")
df <- as.data.frame(GO_results)
write.csv(df,"cl2_GO.csv")

# df_cl1 <- df
# df_cl2 <- df
# 
# subset_df <- df_cl1[df_cl1$ONTOLOGY == "MF", ]
# nrow(subset_df)


KEGG_results <- enrichMKEGG(genes_to_t, organism='hsa', minGSSize=1)
df_KEGG <- as.data.frame(KEGG_results@result)
write.csv(df_KEGG,"cl2_KEGG.csv")


gene_list <- read.csv("KIRP_genes.csv")
genes_to_t <- rownames(gene_list)
GO_results <- enrichGO(gene = genes_to_t, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID", ont = "ALL")
df <- as.data.frame(GO_results)
write.csv(df,"KIRP_GO.csv")


KEGG_results <- enrichMKEGG(genes_to_t, organism='hsa', minGSSize=1)
df_KEGG <- as.data.frame(KEGG_results@result)
write.csv(df_KEGG,"KIRP_KEGG.csv")
