# install necessary package if you have not done so
# install.packages("ggplot2")

# load the library
library(ggplot2)
library(tidyverse)
library(dplyr)
library(matrixStats)
library(cluster)
library(FactoMineR)
library(factoextra)

# Read the data
em_see <- read.csv("em_see.csv")[,-1]
nprobes <- nrow(em_see)
num_e0_probes <- sum(em_see$e0)
em_see$new <- em_see$e0

new <- em_see %>%
  select(master_list, new) %>%
  mutate(new = if_else(new == 0, NA_real_, new)) %>%
  na.omit() %>%
  select(-new) %>%
  setNames("Name")

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
colnames(leftrow) <- "RNAseq_ID"
data <- cbind(leftrow, data)

# Remove rows where all values are NaN in trainingData portion
data <- data[!apply(is.na(data[,(ncol(leftrow)+1):ncol(data)]), 1, all), ]

in_df2 <- data[, 1] %in% em_see$master_list
data <- data[in_df2, ]

# Remove rows where all values are NaN in trainingData portion
filtered_data <- data[data$RNAseq_ID %in% new$Name, ]
filtered_data <- filtered_data[,-1]

notes3 <- read.csv("notes3.csv")
filtered_data <- filtered_data[, notes3$RNA_seq_filename]

t_data <- t(filtered_data)

# Let's scale the data
data_scaled <- scale(t_data)

# Perform PCA
res.pca <- PCA(data_scaled, graph = FALSE)

# # Let's create labels for the groups
# groups <- c(rep("Experimental", 268), rep("Control", 23))
# 
# # Visualize the PCA
# fviz_pca_ind(res.pca,
#              geom.ind = "point",  # show points only (nbut not "text")
#              col.ind = groups, # color by groups
#              palette = c("#00AFBB", "#E7B800"),
#              addEllipses = TRUE, # Concentration ellipses
#              legend.title = "Groups"
# )

# Let's create labels for the groups
groups <- c(rep("Exp1", 50), rep("Exp2", 66), rep("Exp3", 152), rep("Control", 23))

# Visualize the PCA
fviz_pca_ind(res.pca,
             geom.ind = "point",  # show points only (not "text")
             col.ind = groups, # color by groups
             palette = c("#00AFBB", "#E7B800", "#FC4E07", "#00FF00"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups"
)





