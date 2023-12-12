library(dplyr)
library(tidyverse)

bval <- read.csv("bval.csv")

# Calculate the row means
row_means <- rowMeans(bval[,-c(1,2)])

# Create a new data frame with the row means
row_means_df <- data.frame(rowname = bval[, 1], mean = row_means)

# Set the first column as row names in the new data frame
rownames(row_means_df) <- row_means_df$rowname
row_means_df$rowname <- NULL

row_means_df$IlmnID = rownames(row_means_df)

headless <- read.csv("headless.csv")

# Specify the column names to keep
columns_to_keep <- c("IlmnID", "Name", "AddressA_ID", "AddressB_ID", "Genome_Build",
                     "CHR", "MAPINFO", "Chromosome_36", "UCSC_RefGene_Name", "UCSC_RefGene_Accession")

# Assuming your data frame is called 'df', you can use the subset() function to keep only the desired columns
small_headless <- subset(headless, select = columns_to_keep)

# Use the merge function to perform the merge
merged_df <- merge(row_means_df, small_headless, by = "IlmnID", all.x = TRUE)

write.csv(merged_df,"t21.csv")


t26B <- read_csv("t26.csv")
t26B[is.na(t26B)] <- ""
rownames(t26B) <- t26B[,1]
t26B <- t26B[,-c(1,2)]
write.csv(t26B,"t26.csv")
