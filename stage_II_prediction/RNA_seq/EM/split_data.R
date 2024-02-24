# Required library
library(dplyr)
library(matrixStats)

# Load your data
data <- read.csv("tpm_mod_RNA_seq.csv", header = TRUE)
leftrow <- data[,1:3]
data <- data[,-c(1:3)]

filenames <- read.csv("filenames.csv")
colnames(data) <- filenames$names

# Load the file that contains the column names to be retained in the data
filename_tainning <- read.csv("filename_trainning.csv")

# Retrieve the names to be retained from the 'RNA_seq_filename' column
cols_to_keep <- filename_tainning$RNA_seq_filename

# Include the first column of the original data
cols_to_keep <- c(names(data)[1], cols_to_keep)

# Create new data frame with selected columns
trainingData <- data %>% select(all_of(cols_to_keep))


# Log transform the data frame, adding a small constant to avoid taking the log of zero
trainingData <- log2(trainingData + 0.0001)

# Function to calculate z-scores
z_score <- function(x) {
  (x - mean(x)) / sd(x)
}

# Apply z-score function row-wise
trainingData <- t(apply(trainingData, 1, z_score))

# If you'd like to convert your matrix back into a data frame
trainingData <- as.data.frame(trainingData)

# Combine the leftrow and trainingData before removing NaN rows
combinedData <- cbind(leftrow, trainingData)

# Remove rows where all values are NaN in trainingData portion
combinedData <- combinedData[!apply(is.na(combinedData[,(ncol(leftrow)+1):ncol(combinedData)]), 1, all), ]

# Write the data frame to a CSV file
write.csv(combinedData, "trainingData.csv", row.names = TRUE)

