library(tidyverse)

# Read the em_see file and other data files
em_see <- read.csv("em_see.csv")
bval_all <- read_csv("beta_val_all.csv")
mval_all <- read_csv("mval_all.csv")
filenames <- read.csv("filenames.csv")
ref <- read.csv("ref4.csv")

em_see <- em_see[,-1]
em_see[em_see == 0] <- NA
new <- na.omit(em_see)
new <- as.data.frame(new[,1])
colnames(new) <- "Name"

beta_ex <- merge(x = new, y = bval_all, by.x = "Name", by.y = "...1")
mx <- merge(x = new, y = mval_all, by.x = "Name", by.y = "...1")

# pre-processing the files to remove any unwanted column taken up on
# the sql inner join and re-enter corrected files-names for column names
# then transposing the matrix
rownames(mx) = mx[,1]
mx <- mx[,-1]
filenames <- read.csv("filenames.csv")
colnames(mx) <- t(filenames)

# Set seed and perform kmeans clustering on mval_all
set.seed(20)
clusters <- as.data.frame(filenames)
clusters$cluster <- kmeans(t(mx),3,iter.max=10000)$cluster

# Merge ref with clusters
notes <- merge(x = ref, y = clusters, by.x = "File_Name", by.y = "names")
notes1 <- notes[,-c(4, 5, 9, 10, 11, 13, 14, 15)]

# Create tables for each cluster and normal samples
cluster1 <- filter(notes1, cluster == 1, Sample_Type == "Primary Tumor")
cluster2 <- filter(notes1, cluster == 2, Sample_Type == "Primary Tumor")
cluster3 <- filter(notes1, cluster == 3, Sample_Type == "Primary Tumor")
normal <- filter(notes1, Sample_Type == "Solid Tissue Normal")

################### stats.R ######################

print("pre-processing and loading up bval_all...")

# Use bval_all for statistics instead of bx_current
beta_current <- bval_all

temp <- beta_current[,1]
temp <- as.list(t(temp))
beta_current <- beta_current[,-c(1,2)]
rownames(beta_current) <- c(temp)
colnames(beta_current) <- t(filenames)

bval_t <- t(beta_current)
temp <- as.data.frame(rownames(bval_t))
colnames(temp) <- "file_name"
bval_t <- cbind(temp, bval_t)

### Separating information into data frames for statistical info (pre-processing) ###
print("building a data frame for normal + each cluster")

# Merging and transposing for each data frame
# Optimized code for creating normalB data frame
normalB <- merge(x = normal["File_Name"], y = bval_t, by.x = "File_Name", by.y = "file_name", all = FALSE)
rownames(normalB) <- normalB$File_Name
normalB <- t(normalB[,-1])  # Remove the 'File_Name' column and transpose

# Optimized code for creating cluster1B data frame
cluster1B <- merge(x = cluster1["File_Name"], y = bval_t, by.x = "File_Name", by.y = "file_name", all = FALSE)
rownames(cluster1B) <- cluster1B$File_Name
cluster1B <- t(cluster1B[,-1])  # Remove the 'File_Name' column and transpose

# Optimized code for creating cluster2B data frame
cluster2B <- merge(x = cluster2["File_Name"], y = bval_t, by.x = "File_Name", by.y = "file_name", all = FALSE)
rownames(cluster2B) <- cluster2B$File_Name
cluster2B <- t(cluster2B[,-1])  # Remove the 'File_Name' column and transpose

# Optimized code for creating cluster3B data frame
cluster3B <- merge(x = cluster3["File_Name"], y = bval_t, by.x = "File_Name", by.y = "file_name", all = FALSE)
rownames(cluster3B) <- cluster3B$File_Name
cluster3B <- t(cluster3B[,-1])  # Remove the 'File_Name' column and transpose


# Function to calculate row variances
rowVars <- function(x, na.rm=F) {
  rowSums((x - rowMeans(x, na.rm=na.rm))^2, na.rm=na.rm) / (ncol(x) - 1)
}

### Getting statistics for norm + clusters ###
print("getting stats for norm + clusters...")


### Continuing with statistical analysis ###

# Calculating averages for each group
summary <- as.data.frame(rowMeans(normalB))
colnames(summary) <- "norm_avg"
summary$c1_avg <- rowMeans(cluster1B)
summary$c2_avg <- rowMeans(cluster2B)
summary$c3_avg <- rowMeans(cluster3B)

# Performing Kolmogorov–Smirnov tests to get p-values
summary$pval_norm_vs_c1 <- sapply(1:nrow(normalB), function(i) ks.test(as.vector(normalB[i,]), as.vector(cluster1B[i,]))$p.value)
summary$pval_norm_vs_c2 <- sapply(1:nrow(normalB), function(i) ks.test(as.vector(normalB[i,]), as.vector(cluster2B[i,]))$p.value)
summary$pval_norm_vs_c3 <- sapply(1:nrow(normalB), function(i) ks.test(as.vector(normalB[i,]), as.vector(cluster3B[i,]))$p.value)

# Adjusting p-values using Benjamini & Hochberg method
summary$padj_norm_vs_c1 <- p.adjust(summary$pval_norm_vs_c1, method = "BH")
summary$padj_norm_vs_c2 <- p.adjust(summary$pval_norm_vs_c2, method = "BH")
summary$padj_norm_vs_c3 <- p.adjust(summary$pval_norm_vs_c3, method = "BH")

# Writing the summary to a CSV file
write.csv(summary, "bval_all_stat_sum.csv")

# # Creating a gene list from the manifest file
# manifest <- read.csv("headless.csv")
# 
# # Merging with the gene list and separating gene names
# temp <- merge(x = as.data.frame(bval_all[,1]), y = manifest, by.x="V1", by.y="IlmnID", all = FALSE)
# temp <- cbind(temp$V1, temp$UCSC_RefGene_Name)
# colnames(temp) <- c("probes","genes")
# temp <- data.frame(temp)
# 
# # Separating genes and creating a unique gene list
# gene <- temp %>%
#   separate_rows(genes, sep=";")
# genes <- unique(gene)
# 
# # Writing the gene list to a CSV file
# # write.csv(genes, "genes_all.csv")
# 
# print("Analysis complete.")
