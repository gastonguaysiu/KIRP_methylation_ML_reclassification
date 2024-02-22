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

# Prepare the 'new' dataframe
new <- em_see %>%
  select(master_list, new) %>%
  mutate(new = if_else(new == 0, NA_real_, new)) %>%
  na.omit() %>%
  select(-new) %>%
  setNames("Name")

# Load your data
data <- read_csv("rpm_RNA_seq2.csv")
filenames <- read.csv("filenames.csv")

# Retrieve the names to be retained from the 'miRNA_filename' column
# Include the first column of the original data
# Create new data frame with selected columns
cols_to_keep <- filenames$miRNA_filename
cols_to_keep <- c(names(data)[1], cols_to_keep)
data <- data %>% select(all_of(cols_to_keep))

names <- data[,1]
data <- data[,-1]

# Log transform the data frame, adding a small constant to avoid taking the log of zero
data <- log2(data + 0.0001)

# Function to calculate z-scores
z_score <- function(x) {
  (x - mean(x)) / sd(x)
}

# Apply z-score function row-wise
data <- t(apply(data, 1, z_score))

# If you'd like to convert your matrix back into a data frame
data <- as.data.frame(data)

# Combine the leftrow and trainingData before removing NaN rows
data <- cbind(names, data)

# Remove rows where all values are NaN in trainingData portion
data <- data[!apply(is.na(data[,(ncol(names)+1):ncol(data)]), 1, all), ]
rownames(data) <- data[,1]
filtered_data <- data[data$miRNA_ID %in% new$Name, ]
filtered_data <- filtered_data[,-1]

notes3 <- read.csv("notes3.csv")
filtered_data <- filtered_data[, notes3$miRNA_filename]


# Transpose data (make rows into columns and vice versa)
data_transposed <- t(filtered_data)

# Define groups
group <- c(rep('exp1', 92), rep('exp2', 144 - 92), rep('exp3', 236 - 144), rep('control', 291 - 236))

# Perform PCA
pca_res <- PCA(data_transposed, scale.unit = TRUE, graph = FALSE)

# Create a data frame for PCA results and group labels
df <- data.frame(pca_res$ind$coord, group = factor(group))

# Calculate percentage of variance explained by each principal component
explained_var <- round(pca_res$eig[, 2], 1)

# Make a PCA plot with confidence ellipses
ggplot(df, aes(x = Dim.1, y = Dim.2, color = group)) +
  geom_point(alpha = 0.6, size = 2.5) +
  stat_ellipse(level = 0.95) + # Add 95% confidence ellipses
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07", "#0073C2")) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(color = "Group", 
       x = paste("PC1 (", explained_var[1], "%)", sep = ""),
       y = paste("PC2 (", explained_var[2], "%)", sep = ""),
       title = "PCA Plot")
