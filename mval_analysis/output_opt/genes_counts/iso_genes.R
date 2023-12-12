library("tidyverse")
summary <- read.csv("mval_sum.csv")
rownames(summary) <- summary[,1]
summary <- summary[,-1]

for (i in 1:3) {
  cluster_name <- paste0("c", i)
  probes_data <- data.frame(diff = summary[[paste0(cluster_name, "_avg")]] - summary$norm_avg,
                            padj = summary[[paste0("padj_norm_vs_", cluster_name)]])
  row.names(probes_data) <- row.names(summary) # Add row names from summary data frame
  assign(paste0("probes_cl", i), probes_data)
}

###################   sig_probe   ######################

std_diff <- max(sd(probes_cl1$diff), sd(probes_cl2$diff), sd(probes_cl3$diff))

get_sig_probes <- function(cluster, threshold = 0.01) {
  cluster[(cluster$diff > std_diff | cluster$diff < -std_diff) & cluster$padj < threshold,]
}

# Get significant probes
cl1B <- get_sig_probes(probes_cl1)
cl2B <- get_sig_probes(probes_cl2)
cl3B <- get_sig_probes(probes_cl3)

cl1B$probes <- rownames(cl1B)
cl2B$probes <- rownames(cl2B)
cl3B$probes <- rownames(cl3B)

# Read manifest file
manifest <- read_csv("headless.csv")


process_cluster <- function(cluster_data, manifest) {
  cl_temp <- merge(x = cluster_data, y = manifest, by.x = "probes", by.y = "IlmnID", all = FALSE)
  cl_genes <- cl_temp %>%
    separate_rows(UCSC_RefGene_Name, sep = ";")
  unique_cl_genes <- cl_genes %>%
    distinct(UCSC_RefGene_Name, .keep_all = TRUE)
  unique_cl_genes <- unique_cl_genes[,c("probes", "diff", "padj", "UCSC_RefGene_Name")]
  colnames(unique_cl_genes) <- c("probes", "diff", "padj", "gene")
  return(unique_cl_genes)
}

######  sumarize into data frames  ######

cl1_probes_genes <- process_cluster(cl1B, manifest)
cl2_probes_genes <- process_cluster(cl2B, manifest)
cl3_probes_genes <- process_cluster(cl3B, manifest)

cl1_genes <- unique(as.data.frame(cl1_probes_genes$gene))
cl2_genes <- unique(as.data.frame(cl2_probes_genes$gene))
cl3_genes <- unique(as.data.frame(cl3_probes_genes$gene))
colnames(cl1_genes) <- colnames(cl2_genes) <- colnames(cl3_genes) <- "gene"

######################    building simp heatmap data    ######################

sig_probes <- as.data.frame(unique(c(cl1B$probes, cl2B$probes, cl3B$probes)))
colnames(sig_probes) <- "probes"

temp <- summary
temp$probes <- rownames(summary)
temp2 <- merge(x=sig_probes, y=temp, by ="probes", all = FALSE)

simpheat <- temp2[,1:5]
rownames(simpheat) <- simpheat[,1]
simpheat <- simpheat[,-1]

write.csv(sig_probes,"sig_probes.csv")
write.csv(simpheat,"simp_heat.csv")

######  print out the csv  ######

# write.csv(cl1_genes,"cl1_genes.csv", row.names = FALSE)
# write.csv(cl2_genes,"cl2_genes.csv", row.names = FALSE)
# write.csv(cl3_genes,"cl3_genes.csv", row.names = FALSE)

kirp_genes <- read.csv("KIRP_genes.csv")

cl1_kirp <- inner_join(kirp_genes, cl1_genes, by = c("Symbol" = "gene"))
cl2_kirp <- inner_join(kirp_genes, cl2_genes, by = c("Symbol" = "gene"))
cl3_kirp <- inner_join(kirp_genes, cl3_genes, by = c("Symbol" = "gene"))

cl1_kirp <- as.data.frame(cl1_kirp$Symbol)
cl2_kirp <- as.data.frame(cl2_kirp$Symbol)
cl3_kirp <- as.data.frame(cl3_kirp$Symbol)



# save.image("my_workspace.RData")
load("my_workspace.RData")
