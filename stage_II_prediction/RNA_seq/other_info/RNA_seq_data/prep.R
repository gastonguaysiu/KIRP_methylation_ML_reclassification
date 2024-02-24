install.packages("readr")
install.packages("dplyr")
install.packages("purrr")
install.packages("matrixStats")

library(readr)
library(dplyr)
library(purrr)
library(matrixStats)

# Create new directory if it doesn't exist
new_dir <- "new_directory"
if (!dir.exists(new_dir)) {
  dir.create(new_dir)
}

# Get list of all .tsv files in the working directory
tsv_files <- list.files(pattern = "\\.tsv$")

# Function to load .tsv file, modify it and save to new directory
process_tsv <- function(file) {
  # Load the file, skipping the first row, then set the header names with second row
  df <- read_tsv(file, skip = 1, col_names = TRUE)

  # Remove the first 5 rows
  df <- df[-c(1:4), ]

  # Save the modified data frame to a new .tsv file in the new directory
  write_tsv(df, file.path(new_dir, file))
}

# Apply the function to each .tsv file
purrr::walk(tsv_files, process_tsv)


############################ RNA-seq ######################################


# Get list of all .tsv files in the working directory
tsv_files <- list.files(pattern = "\\.tsv$")

# Initialize an empty data frame to store the unique data
unique_df <- tibble()

# Loop over the tsv files
for (i in seq_along(tsv_files)) {
  
  # Read in the tsv file
  df <- read_tsv(tsv_files[i])
  
  # If this is the first file, create a new data frame with just the first three columns
  if (i == 1) {
    unique_df <- df %>%
      select(gene_id, gene_name, gene_type) %>%
      distinct()
  }
  
  # Add the tpm_unstranded column to the unique_df data frame, named after the current file
  unique_df <- unique_df %>%
    left_join(df %>%
                select(gene_id, tpm_unstranded) %>%
                rename(!!paste0(tools::file_path_sans_ext(tsv_files[i]), ".tsv") := tpm_unstranded),
              by = "gene_id")
  
  # Clear the df to free up memory
  rm(df)
  gc()  # Force garbage collection
}

# Replace all NA values with the median of the respective row (only consider numeric columns)
unique_df <- unique_df %>% 
  mutate(across(where(is.numeric), ~ifelse(is.na(.), rowMedians(as.matrix(unique_df[, sapply(unique_df, is.numeric)]), na.rm = TRUE), .), .names = "{.col}"))


write.csv(unique_df,"tpm_mod_RNA_seq.csv")


######################## miRNA-seq ###################################


# Get list of all .tsv files in the working directory
tsv_files <- list.files(pattern = "\\.tsv$")

# Initialize an empty data frame to store the unique data
unique_df <- tibble()

# Loop over the tsv files
for (i in seq_along(tsv_files)) {
  
  # Read in the tsv file
  df <- read_tsv(tsv_files[i])
  
  # If this is the first file, create a new data frame with just the miRNA_ID column
  if (i == 1) {
    unique_df <- df %>%
      select(miRNA_ID) %>%
      distinct()
  }
  
  # Add the reads_per_million_miRNA_mapped column to the unique_df data frame, named after the current file
  unique_df <- unique_df %>%
    left_join(df %>%
                select(miRNA_ID, reads_per_million_miRNA_mapped) %>%
                rename(!!paste0(tools::file_path_sans_ext(tsv_files[i]), ".tsv") := reads_per_million_miRNA_mapped),
              by = "miRNA_ID")
  
  # Clear the df to free up memory
  rm(df)
  gc()  # Force garbage collection
}

# Replace all NA values with the median of the respective row (only consider numeric columns)
unique_df <- unique_df %>% 
  mutate(across(where(is.numeric), ~ifelse(is.na(.), rowMedians(as.matrix(unique_df[, sapply(unique_df, is.numeric)]), na.rm = TRUE), .), .names = "{.col}"))


write.csv(unique_df,"rpm_mod_RNA_seq.csv")
