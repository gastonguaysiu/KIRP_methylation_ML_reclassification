# Load necessary libraries
library(tidyverse)
library(dplyr)
library(matrixStats)
library(cluster)

# Read the data
filenames <- read.csv("filename_trainning.csv")
miRNA <- read_csv("trainingData.csv")
em_see <- read.csv("em_see.csv")[,-1]

# Compute some quantities
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

# Merge the dataframes
mx <- merge(x = new, y = miRNA, by.x = "Name", by.y = "...1") %>%
  select(-1)

# Set the seed for reproducibility
set.seed(42)

# Perform k-means clustering
kmeans_result <- kmeans(t(mx), 3, iter.max = 10000)

# Save cluster assignments
clusters <- tibble(sample = colnames(mx), cluster = kmeans_result$cluster)

# Save the coordinates of the centroids
centroids <- kmeans_result$centers


# Create a function to find the nearest centroid for a new point
assign_to_cluster <- function(point, centroids) {
  distances <- sqrt(rowSums((t(t(centroids) - point))^2))
  return(which.min(distances))
}

# Apply the function to each new data point
cluster_assignments <- apply(new_data, 1, assign_to_cluster, centroids = loaded_centroids)





# Load your data
data <- read.csv("rpm_RNA_seq.csv", row.names = 1, header = TRUE)
filenames <- read.csv("filenames.csv")
colnames(data) <- filenames$filename

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

# Remove rows where all values are NaN
data <- data[!apply(is.na(data), 1, all), ]

# Write the data frame to a CSV file
write.csv(data, "temp.csv", row.names = TRUE)


filenames_trainning <- read.csv("filename_trainning.csv")
miRNA <- read_csv("temp.csv")

em_see <- read.csv("em_see.csv")[,-1]
nprobes <- nrow(em_see)
num_e0_probes <- sum(em_see$e0)

nprobes <- nrow(em_see)
num_e0_probes <- sum(em_see$e0)
em_see$new <- em_see$e0

new <- em_see %>%
  select(master_list, new) %>%
  mutate(new = if_else(new == 0, NA_real_, new)) %>%
  na.omit() %>%
  select(-new) %>%
  setNames("Name")

mx <- merge(x = new, y = miRNA, by.x = "Name", by.y = "...1") %>%
  select(-1)

set.seed(42)
cluster_assignments <- apply(data, 1, assign_to_cluster, centroids = centroids)









ref <- read.csv("gdc_join_clinical.csv")
notes <- merge(x = ref, y = clusters, by.x = "miRNA_filename", by.y = "sample") %>%
  select(-c(5,6,15:28))

notes2 <- anti_join(notes, filenames_trainning, by = "miRNA_filename")

calculate_scores <- function(df, cluster_id) {
  df <- df %>%
    filter(cluster == cluster_id, Sample_Type == "Primary Tumor") %>%
    mutate(live_score = as.numeric(days_to_last_follow_up) / 624,
           death_score = as.numeric(days_to_death) / 761)
}


mxcls <- lapply(1:3, calculate_scores, df = notes2)

nx <- data.frame(samples = vapply(mxcls, nrow, numeric(1)),
                 dead = vapply(mxcls, function(df) sum(df$vital_status == "Dead"), numeric(1)),
                 live_weight = vapply(mxcls, function(df) sum(df$live_score, na.rm = TRUE), numeric(1)),
                 dead_weight = 0.01)

total_dead <- sum(nx$dead)

nx <- nx %>%
  rownames_to_column("cluster") %>%
  mutate(dead = if_else(dead == 0, 0.01, dead),
         dead_weight = if_else(dead > 0, if_else(vapply(mxcls, function(df) sum(df$death_score, na.rm = TRUE), numeric(1)) == 0, 0.01, vapply(mxcls, function(df) sum(df$death_score, na.rm = TRUE), numeric(1))), dead_weight),
         ratio = dead / samples,
         precision = dead / (dead + live_weight),
         recall = dead / total_dead,
         final_score = if_else((precision + recall) == 0, 0, 2 * (precision * recall) / (precision + recall)),
         rscore = ratio * dead_weight)

nxa <- nx %>% arrange(final_score)
nxb <- nx %>% arrange(rscore)
write.csv(nxa, "na_test.csv")
write.csv(nxb, "nb_test.csv")


mxcls <- lapply(1:3, calculate_scores, df = notes)

nx <- data.frame(samples = vapply(mxcls, nrow, numeric(1)),
                 dead = vapply(mxcls, function(df) sum(df$vital_status == "Dead"), numeric(1)),
                 live_weight = vapply(mxcls, function(df) sum(df$live_score, na.rm = TRUE), numeric(1)),
                 dead_weight = 0.01)

total_dead <- sum(nx$dead)

nx <- nx %>%
  rownames_to_column("cluster") %>%
  mutate(dead = if_else(dead == 0, 0.01, dead),
         dead_weight = if_else(dead > 0, if_else(vapply(mxcls, function(df) sum(df$death_score, na.rm = TRUE), numeric(1)) == 0, 0.01, vapply(mxcls, function(df) sum(df$death_score, na.rm = TRUE), numeric(1))), dead_weight),
         ratio = dead / samples,
         precision = dead / (dead + live_weight),
         recall = dead / total_dead,
         final_score = if_else((precision + recall) == 0, 0, 2 * (precision * recall) / (precision + recall)),
         rscore = ratio * dead_weight)

nxa <- nx %>% arrange(final_score)
nxb <- nx %>% arrange(rscore)
write.csv(nxa, "na_all.csv")
write.csv(nxb, "nb_all.csv")
