# Loading necessary libraries
library(dplyr)
library(readr)

# Setting the seed for reproducibility
set.seed(42)

# Reading in the data
gdc_data <- read_csv("gdc_join_clinical.csv")

# Calculate the number of rows for training set (90% of total)
train_size <- floor(0.9 * nrow(gdc_data))

# Randomly sample rows for the training set
train_indices <- sample(seq_len(nrow(gdc_data)), size = train_size)

# Create the training set
training_data <- gdc_data[train_indices, ]

# Writing the training and testing sets to CSV files
write_csv(training_data, "gdc_train.csv")
