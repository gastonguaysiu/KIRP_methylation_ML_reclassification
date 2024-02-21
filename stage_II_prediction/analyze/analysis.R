# Load necessary libraries
library(tidyverse)
library(dplyr)
library(matrixStats)
library(cluster)

wssplot <- function(f1, nc=5, seed=1234){
  wss <- (nrow(f1)-1)*sum(apply(f1,2,var))
  for (i in 2:nc){
    set.seed(seed)
    wss[i] <- sum(kmeans(f1, centers=i)$withinss)}
  plot(1:nc, wss, type="b", xlab="Number of Clusters",
       ylab="Within groups sum of squares")
  wss
}

# Read the data
RNA_seq <- read_csv("trainingData.csv")
em_see <- read.csv("em_see.csv")[,-1]
filename_trainning <- read.csv("filename_trainning.csv")

# Compute some quantities
nprobes <- nrow(em_see)
num_e0_probes <- sum(em_see$e0)
em_see$new <- em_see$e0

# Prepare the 'new' dataframe
new <- em_see %>%
  select(master_list, new) %>%
  mutate(new = if_else(new == 0, NA_real_, new)) %>%
  na.omit() %>%
  select(-new) %>%
  setNames("Name")

# Merge the dataframes
mx <- merge(x = new, y = RNA_seq, by.x = "Name", by.y = "...1") %>%
  select(-1)

# Set the seed for reproducibility
set.seed(42)
kmeans_result <- kmeans(t(mx), 3, iter.max = 10000)
clusters <- tibble(sample = colnames(mx), cluster = kmeans_result$cluster)
centroids <- kmeans_result$centers


# Load your data
data <- read.csv("tpm_RNA_seq2.csv", header = TRUE)
leftrow <- as.data.frame(data[,1])
data <- data[,-1]

filenames <- read.csv("filenames.csv")
colnames(data) <- filenames$names

filenames <- read.csv("filenames2.csv")

# Retrieve the names to be retained from the 'RNA_seq_filename' column
# Include the first column of the original data
# Create new data frame with selected columns
cols_to_keep <- filenames$RNA_seq_filename
cols_to_keep <- c(names(data)[1], cols_to_keep)
data <- data %>% select(all_of(cols_to_keep))

# Log transform the data frame, adding a small constant to avoid taking the log of zero
data <- log2(data + 0.0001)

# calculate row means
# replace NA values with row means
row_means <- rowMeans(data, na.rm = TRUE)
data[is.na(data)] <- row_means[row(data)[is.na(data)]]

# Function to calculate z-scores
z_score <- function(x) {
  (x - mean(x)) / sd(x)
}

# Apply z-score function row-wise
data <- t(apply(data, 1, z_score))

# If you'd like to convert your matrix back into a data frame
data <- as.data.frame(data)

# Combine the leftrow and trainingData before removing NaN rows
data <- cbind(leftrow, data)

# Remove rows where all values are NaN in trainingData portion
data <- data[!apply(is.na(data[,(ncol(leftrow)+1):ncol(data)]), 1, all), ]

in_df2 <- data[, 1] %in% em_see$master_list
data <- data[in_df2, ]

# Merge the dataframes
mx <- merge(x = new, y = data, by.x = "Name", by.y = "data[, 1]") %>%
  select(-1)

# Create a function to find the nearest centroid for a new point
assign_to_cluster <- function(point, centroids) {
  distances <- sqrt(rowSums((t(t(centroids) - point))^2))
  return(which.min(distances))
}

# Apply the function to each new data point
cluster_assignments <- apply(mx, 2, assign_to_cluster, centroids = centroids)
new_clusters <- tibble(sample = colnames(mx), cluster = cluster_assignments)

ref <- read.csv("gdc_join_clinical.csv")
notes <- merge(x = ref, y = new_clusters, by.x = "RNA_seq_filename", by.y = "sample") %>%
  select(-c(5,6,15:28))

notes2 <- anti_join(notes, filename_trainning, by = "RNA_seq_filename")




pdf("elbow1.pdf")
wssplot(mx, nc=10)
dev.off()


# Separating information into data frames to later get statistical info (pre-processing)
print("building a data frame for normal + each cluster")

normal <- filter(notes, Sample_Type == "Solid Tissue Normal")
cluster1 <- filter(notes, cluster == 1, Sample_Type == "Primary Tumor")
cluster2 <- filter(notes, cluster == 2, Sample_Type == "Primary Tumor")
cluster3 <- filter(notes, cluster == 3, Sample_Type == "Primary Tumor")
###################   stats.R   ######################


print("pre-processing and loading up mx_current...")

probe_list <- as.data.frame(data[,1])
mx_t <- as.data.frame(lapply(data.frame(t(data), stringsAsFactors = FALSE), as.numeric))
mx_t <- mx_t[-1,]
colnames(mx_t) <- t(probe_list)
temp <- as.data.frame(colnames(data[,-1]))
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

normalB <- build_df(data.frame(sample = normal$RNA_seq_filename))
cluster1B <- build_df(data.frame(sample = cluster1$RNA_seq_filename))
cluster2B <- build_df(data.frame(sample = cluster2$RNA_seq_filename))
cluster3B <- build_df(data.frame(sample = cluster3$RNA_seq_filename))

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

write.csv(summary, "stat_RNAseq_sum.csv")

###################   sig_probe   ######################

for (i in 1:3) {
  cluster_name <- paste0("c", i)
  probes_data <- data.frame(diff = summary[[paste0(cluster_name, "_avg")]] - summary$norm_avg,
                            padj = summary[[paste0("padj_norm_vs_", cluster_name)]])
  row.names(probes_data) <- row.names(summary) # Add row names from summary data frame
  assign(paste0("probes_cl", i), probes_data)
}

std_diff <- max(sd(probes_cl1$diff),sd(probes_cl2$diff), sd(probes_cl3$diff))

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

######################    building simp heatmap data    ######################

sig_probes <- unique(c(cl1B$probes, cl2B$probes, cl3B$probes))
sig_probes <- as.data.frame(sig_probes)
colnames(sig_probes) <- "probes"

get_cluster_means <- function(cluster_data) {
  rowMeans(cluster_data)
}

simp_heat <- data.frame(normal_avg = get_cluster_means(normalB),
                        cl1_avg = get_cluster_means(cluster1B),
                        cl2_avg = get_cluster_means(cluster2B),
                        cl3_avg = get_cluster_means(cluster3B))

simp_heat$temp <- rownames(simp_heat)
simp_heat2 <- merge(x=sig_probes, y=simp_heat, by.x="probes", by.y="temp", all = FALSE)
rownames(simp_heat2) <- simp_heat2[,1]
simp_heat2 <- simp_heat2[,-1]

notes3 <- notes[order(notes$cluster),]

write.csv(notes,"notes.csv")
write.csv(notes2, "notes2.csv")
write.csv(notes3, "notes3.csv")
write.csv(simp_heat,"simp_heat.csv")
write.csv(simp_heat2,"simp_heat2.csv")


save.image("my_workspace.RData")
# load("my_workspace.RData")

