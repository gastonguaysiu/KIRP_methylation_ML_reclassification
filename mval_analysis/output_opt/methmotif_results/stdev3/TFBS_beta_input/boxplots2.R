# Load necessary libraries
library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
library(stats) # For ANOVA and TukeyHSD

# Load bval_all_stat_sum.csv
bval_all_stat_sum <- read_csv("bval_all_stat_sum.csv")

# Load TF_input_regions_mani_v3.csv
TF_input_regions_mani_v3 <- read_csv("TF_input_regions_mani_v3.csv")

# Extract the last two columns of TF_input_regions_mani_v3
TF_last_two_columns <- select(TF_input_regions_mani_v3, 2, 35, 36)

# Perform an inner join using the "Name" column
joined_df <- inner_join(bval_all_stat_sum, TF_last_two_columns, by = "Name")

# Split joined_df into a list of data frames, each corresponding to a unique TF value
list_of_dfs <- split(joined_df, joined_df$TF)

# Initialize data frames for ANOVA p-values
anova_results <- data.frame(ConditionPair = character(), TF = character(), PValue = numeric(), stringsAsFactors = FALSE)

# Define all possible pairs of conditions
condition_pairs <- combn(c("vivo", "HEK293", "c1_avg", "c2_avg", "c3_avg"), 2, simplify = FALSE)

# Loop through each TF-specific data frame
for(tf_name in names(list_of_dfs)) {
  current_df <- list_of_dfs[[tf_name]]
  
  # Check if the current data frame has more than 10 rows
  if(nrow(current_df) > 15) {
    current_df <- current_df %>%
      rename(vivo = norm_avg, HEK293 = beta_val) %>%
      mutate(HEK293 = HEK293 / 100)
    
    # Convert data from wide to long format for analysis
    df_long <- pivot_longer(current_df, 
                            cols = c(vivo, HEK293, c1_avg, c2_avg, c3_avg),
                            names_to = "Condition", 
                            values_to = "Value")
    
    # ANOVA for each condition pair
    for(pair in condition_pairs) {
      condition_formula <- paste("Value ~ Condition", sep = "")
      anova_model <- aov(as.formula(condition_formula), data = df_long %>% filter(Condition %in% pair))
      p_value <- summary(anova_model)[[1]]$"Pr(>F)"[1]
      anova_results <- rbind(anova_results, data.frame(ConditionPair = paste(pair, collapse = " vs "), TF = tf_name, PValue = p_value))
    }
    
  }
}

# At this point, anova_results data frames contain the ANOVA p-values

# Apply the Benjamini-Hochberg procedure to adjust p-values for multiple comparisons
anova_results <- anova_results %>%
  group_by(TF) %>%
  mutate(FDR = p.adjust(PValue, method = "BH")) %>%
  ungroup()

# Optionally, filter to keep only significant results based on a common threshold (e.g., FDR < 0.05)
significant_results <- filter(anova_results, FDR < 0.05)

# Pivot for PValues
p_values_wide <- anova_results %>%
  select(ConditionPair, TF, PValue) %>%
  pivot_wider(names_from = TF, values_from = PValue)

# Pivot for FDRs
fdr_values_wide <- anova_results %>%
  select(ConditionPair, TF, FDR) %>%
  pivot_wider(names_from = TF, values_from = FDR)

# Rename the first column to "Condition" for clarity if needed
colnames(p_values_wide)[1] <- "Condition"
colnames(fdr_values_wide)[1] <- "Condition"

