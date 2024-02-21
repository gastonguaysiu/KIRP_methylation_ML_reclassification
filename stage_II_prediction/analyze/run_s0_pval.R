library(tidyverse)

wssplot <- function(f1, nc=5, seed=1234){
  wss <- (nrow(f1)-1)*sum(apply(f1,2,var))
  for (i in 2:nc){
    set.seed(seed)
    wss[i] <- sum(kmeans(f1, centers=i)$withinss)}
  plot(1:nc, wss, type="b", xlab="Number of Clusters",
       ylab="Within groups sum of squares")
  wss
}

em_see <- read.csv("em_see.csv")
em_see <- em_see[,-1]
em_see[em_see == 0] <- NA
new <- na.omit(em_see)
new <- as.data.frame(new[,1])
colnames(new) <- "Name"

bval_all <- read_csv("bval_all.csv")
mval_all <- read_csv("mval_all.csv")
probe_list <- as.data.frame(mval_all$"...1")

bx <- merge(x = new, y = bval_all, by.x = "Name", by.y = "...1")
mx <- merge(x = new, y = mval_all, by.x = "Name", by.y = "...1")

rownames(mx) <- mx[,1]
mx <- mx[,-1]

pdf("elbow1.pdf")
wssplot(mx, nc=10)
dev.off()

set.seed(20)
clusters <- as.data.frame(colnames(mx))
clusters$cluster <- kmeans(t(mx), 3, iter.max = 10000)$cluster
colnames(clusters) <- c("sample", "cluster")

ref <- read.csv("n2.csv")
notes <- inner_join(ref, clusters, by = "sample")

normal <- notes %>% filter( Braak_stage == "0")
cluster1 <- notes %>% filter(cluster == 1, Braak_stage != "0")
cluster2 <- notes %>% filter(cluster == 2, Braak_stage != "0")
cluster3 <- notes %>% filter(cluster == 3, Braak_stage != "0")

###################   stats.R   ######################
print("pre-processing and loading up mx_current...")

mx_t <- as.data.frame(lapply(data.frame(t(mval_all), stringsAsFactors = FALSE), as.numeric))
mx_t <- mx_t[-1,]
colnames(mx_t) <- t(probe_list)
temp <- as.data.frame(colnames(mval_all[,-1]))
colnames(temp) <- "sample"
mx_t <- cbind(temp, mx_t)

# Function to process clusters
build_df <- function(cluster_samples) {
  cluster_df <- merge(x = cluster_samples, y = mx_t, by = "sample")
  rownames(cluster_df) <- cluster_df$sample
  return(as.data.frame(t(cluster_df[-1])))
}

# Separating information into data frames to later get statistical info (pre-processing)
print("building a data frame for normal + each cluster")

normalB <- build_df(data.frame(sample = normal$sample))
cluster1B <- build_df(data.frame(sample = cluster1$sample))
cluster2B <- build_df(data.frame(sample = cluster2$sample))
cluster3B <- build_df(data.frame(sample = cluster3$sample))

# Function to calculate row variance
rowVars <- function(x, na.rm = F) {
  rowSums((x - rowMeans(x, na.rm = na.rm))^2, na.rm = na.rm) / (ncol(x) - 1)
}

# Function to perform Kolmogorov–Smirnov test of two samples (cluster vs norm) and update summary
update_summary <- function(summary, norm_data, cluster_data, cluster_name) {
  summary[[paste0(cluster_name, "_avg")]] <- rowMeans(cluster_data)
  summary[[paste0(cluster_name, "_std")]] <- sqrt(rowVars(cluster_data))
  
  pval_col_name <- paste0("pval_norm_vs_", cluster_name)
  summary[[pval_col_name]] <- sapply(1:nrow(norm_data), function(i) {
    ks.test(as.numeric(norm_data[i, ]), as.numeric(cluster_data[i, ]))$p
  })
  
  summary[[paste0("padj_norm_vs_", cluster_name)]] <- p.adjust(summary[[pval_col_name]], method = "BH")
  
  return(summary)
}

# Getting stats for norm + clusters
print("getting stats for norm + clusters...")

summary <- data.frame(norm_avg = rowMeans(normalB))
summary$norm_std <- sqrt(rowVars(normalB))

# Update summary for each cluster
cluster_names <- c("c1", "c2", "c3")
cluster_data_list <- list(cluster1B, cluster2B, cluster3B)

for (i in seq_along(cluster_names)) {
  summary <- update_summary(summary, normalB, cluster_data_list[[i]], cluster_names[i])
}

write.csv(summary, "stat_mval_sum.csv")

# Create some tables to export to sqlite3, where we restrict to only significant probes
# i.e difference in m-val expression of 5% and a p-adjusted of less than 0.01
for (i in 1:3) {
  cluster_name <- paste0("c", i)
  probes_data <- data.frame(diff = summary[[paste0(cluster_name, "_avg")]] / summary$norm_avg, pval = summary[[paste0("pval_norm_vs_", cluster_name)]])
  row.names(probes_data) <- row.names(summary) # Add row names from summary data frame
  assign(paste0("probes_cl", i), probes_data)
}

######################    building simp heatmap data    ######################

get_cluster_means <- function(cluster_data) {
  rowMeans(cluster_data)
}

simp_heat <- data.frame(cl1_avg = get_cluster_means(cluster1B),
                        cl2_avg = get_cluster_means(cluster2B),
                        cl3_avg = get_cluster_means(cluster3B))

###################   sig_probe.sh   ######################

get_sig_probes <- function(cluster, threshold = 0.01) {
  cluster[(cluster$diff > 1.05 | cluster$diff < 0.95) & cluster$pval < threshold,]
}

# Get significant probes
cl1B <- get_sig_probes(probes_cl1)
cl2B <- get_sig_probes(probes_cl2)
cl3B <- get_sig_probes(probes_cl3)

# Create data frames with probe names
cl1b <- as.data.frame(rownames(cl1B), stringsAsFactors = FALSE)
cl2b <- as.data.frame(rownames(cl2B), stringsAsFactors = FALSE)
cl3b <- as.data.frame(rownames(cl3B), stringsAsFactors = FALSE)
colnames(cl1b) <- colnames(cl2b) <- colnames(cl3b) <- c("probes")


# Combine significant probes and remove duplicates
sig_probes <- unique(rbind(cl1b, cl2b, cl3b))
colnames(sig_probes) <- c("probes")

# Read manifest file
manifest <- read_csv("headless.csv")

# Define a function to process the data
process_cluster <- function(cluster_data, manifest) {
  cl_temp <- merge(x = cluster_data, y = manifest, by.x = "probes", by.y = "IlmnID", all = FALSE)
  cl_temp <- as.data.frame(cl_temp$UCSC_RefGene_Name)
  colnames(cl_temp) <- c("genes")
  cl_genes <- cl_temp %>%
    separate_rows(genes, sep = ";")
  unique_cl_genes <- unique(cl_genes)
  return(unique_cl_genes)
}

# Process each cluster using the function
cl1_dm_genes <- process_cluster(cl1b, manifest)
cl2_dm_genes <- process_cluster(cl2b, manifest)
cl3_dm_genes <- process_cluster(cl3b, manifest)

######  make some files for GO enrichment analysis  ######

write.table(cl1_dm_genes, "cl1_dm_genes.csv", sep = ",", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(cl2_dm_genes, "cl2_dm_genes.csv", sep = ",", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(cl3_dm_genes, "cl3_dm_genes.csv", sep = ",", col.names = FALSE, row.names = FALSE, quote = FALSE)

save.image("my_workspace.RData")
load("my_workspace.RData")

