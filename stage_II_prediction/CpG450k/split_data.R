# Required library
library(tidyverse)

# Load your data
data <- read_csv("mval_all.csv")
data <- as.data.frame(data)
rownames(data) <- data[,1]
data <- data[,-1]

filenames <- read.csv("filenames.csv")
colnames(data) <- filenames$names

# Load the file that contains the column names to be retained in the data
filenames_training <- read.csv("filename_trainning.csv")

# Retrieve the names to be retained from the 'miRNA_filename' column
cols_to_keep <- filenames_training$CpG_filename

# Include the first column of the original data
cols_to_keep <- c(names(data)[1], cols_to_keep)

# Create new data frame with selected columns
trainingData <- data %>% select(all_of(cols_to_keep))

# Write the data frame to a CSV file
write.csv(trainingData, "trainingData.csv", row.names = TRUE)


