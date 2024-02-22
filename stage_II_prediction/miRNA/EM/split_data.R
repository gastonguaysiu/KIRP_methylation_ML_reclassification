# Required library
library(dplyr)
library(matrixStats)


# Load your data
data <- read.csv("rpm_RNA_seq.csv", row.names = 1, header = TRUE)
filenames <- read.csv("filenames.csv")
colnames(data) <- filenames$filename

# Load the file that contains the column names to be retained in the data
filenames_training <- read.csv("filename_trainning.csv")

# Retrieve the names to be retained from the 'miRNA_filename' column
cols_to_keep <- filenames_training$miRNA_filename

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

# Remove rows where all values are NaN
trainingData <- trainingData[!apply(is.na(trainingData), 1, all), ]

# Write the data frame to a CSV file
write.csv(trainingData, "trainingData.csv", row.names = TRUE)
