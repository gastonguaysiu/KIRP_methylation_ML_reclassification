# Load necessary libraries
library(dplyr)
library(readr)

# Step 1: Conversion of Coordinates
# Load all_TFBS_input_plus3.bed as a TSV file with an additional column for TF
col_names <- c("hg38_chr", "hg38_start", "hg38_end", "beta_val", "TF")
all_TFBS_hg38 <- read_tsv("all_TFBS_input_plus3.bed", col_names = col_names, col_types = cols())

# Load conversion table
conversion <- read.csv("conversion.csv")

# Adjust hg38_start by subtracting 1 to align correctly with conversion table
# all_TFBS_hg38$hg38_start <- all_TFBS_hg38$hg38_start - 1

# Convert hg38 coordinates to hg19 by matching on hg38_chr and hg38_end
converted_TFBS <- all_TFBS_hg38 %>%
  left_join(conversion, by = c("hg38_chr", "hg38_end")) %>%
  select(hg19_chr, hg19_start, hg19_end, beta_val, TF) # Include TF in the selection

# Remove rows with NA in the hg19_chr column
converted_TFBS <- filter(converted_TFBS, !is.na(hg19_chr))
converted_TFBS <- unique(converted_TFBS)

# Save the converted data to a new CSV file
write.csv(converted_TFBS, "all_TFBS_input_hg19.csv", row.names = FALSE)

# Step 2: Matching with Headless Data
# Load headless.csv
headless <- read.csv("headless.csv")

# Load the converted all_TFBS_input_hg19.csv
all_TFBS_hg19 <- read.csv("all_TFBS_input_hg19.csv")

# Adjust the hg19_chr column by removing the "chr" prefix
all_TFBS_hg19$hg19_chr <- gsub("^chr", "", all_TFBS_hg19$hg19_chr)

# Find matching rows
matched_data <- inner_join(headless, all_TFBS_hg19, by = c("CHR" = "hg19_chr", "MAPINFO" = "hg19_end"))

write.csv(matched_data, 'TF_input_regions_mani_v3.csv', row.names = FALSE)
