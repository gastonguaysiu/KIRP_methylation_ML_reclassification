# Load necessary libraries
library(tidyverse)
library(dplyr)
library(matrixStats)
library(cluster)

# Read the data
CpG <- read_csv("trainingData.csv")
em_see <- read.csv("em_see.csv")[,-1]
filename_trainning <- read.csv("filename_trainning.csv")

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
mx <- merge(x = new, y = CpG, by.x = "Name", by.y = "...1") %>%
  select(-1)

# Set the seed for reproducibility
set.seed(42)
kmeans_result <- kmeans(t(mx), 3, iter.max = 10000)
clusters <- tibble(sample = colnames(mx), cluster = kmeans_result$cluster)
centroids <- kmeans_result$centers


ref <- read.csv("gdc_join_clinical.csv")
notes <- merge(x = ref, y = clusters, by.x = "CpG_filename", by.y = "sample") %>%
  select(-c(5,6,15:28))

live_score_constant <- median(ref$days_to_death, na.rm = TRUE)
death_score_constant <- median(ref$days_to_last_follow_up, na.rm = TRUE)

calculate_scores <- function(df, cluster_id) {
  df <- df %>%
    filter(cluster == cluster_id, Sample_Type == "Primary Tumor") %>%
    mutate(live_score = as.numeric(days_to_last_follow_up) / live_score_constant,
           death_score = as.numeric(days_to_death) / death_score_constant)
}


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
write.csv(nxa, "na_best.csv")
write.csv(nxb, "nb_best.csv")

#######################################33

# Load your data
data <- read_csv("mval_all2.csv")
data2 <- as.data.frame(data)
rownames(data2) <- data2[,1]
data2 <- data2[,-1]

filenames <- read.csv("filenames2.csv")

# Merge the dataframes
mx <- merge(x = new, y = data, by.x = "Name", by.y = "...1") %>%
  select(-1)


# Create a function to find the nearest centroid for a new point
assign_to_cluster <- function(point, centroids) {
  distances <- sqrt(rowSums((t(t(centroids) - point))^2))
  return(which.min(distances))
}

# Apply the function to each new data point
cluster_assignments <- apply(mx, 2, assign_to_cluster, centroids = centroids)
new_clusters <- tibble(sample = colnames(mx), cluster = cluster_assignments)


notes <- merge(x = ref, y = new_clusters, by.x = "CpG_filename", by.y = "sample") %>%
  select(-c(5,6,15:28))

notes2 <- anti_join(notes, filename_trainning, by = "CpG_filename")

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


#############################################

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
         