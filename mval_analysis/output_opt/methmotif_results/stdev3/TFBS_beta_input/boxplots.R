# Load necessary libraries
library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)

# Load bval_all_stat_sum.csv
bval_all_stat_sum <- read_csv("bval_all_stat_sum.csv")

# Load TF_input_regions_mani_v3.csv
TF_input_regions_mani_v3 <- read_csv("TF_input_regions_mani_v3.csv")

# Extract the last two columns of TF_input_regions_mani_v3
TF_last_two_columns <- select(TF_input_regions_mani_v3, 2, 35, 36)

# Perform an inner join using the "Name" column, including all columns from bval_all_stat_sum
# and only the last two columns from TF_input_regions_mani_v3
joined_df <- inner_join(bval_all_stat_sum, TF_last_two_columns, by = "Name")

# # Optional: Save the joined data frame to a new CSV file
# write_csv(joined_data_frame, "joined_data_frame.csv")

# Split joined_df into a list of data frames, each corresponding to a unique TF value
list_of_dfs <- split(joined_df, joined_df$TF)

# Optionally, you might want to perform operations on each data frame in the list
# For example, to simply print out the number of rows for each TF's data frame, you could do:
lapply(list_of_dfs, function(df) nrow(df))

# If you need to access a specific data frame for a TF, you can do so by:
# specific_df <- list_of_dfs[["specific_TF_value"]]
# Replace "specific_TF_value" with the actual TF value you're interested in.

# # For example, to save each data frame to a separate CSV file named after the TF value:
# lapply(names(list_of_dfs), function(tf_name) {
#   write_csv(list_of_dfs[[tf_name]], paste0(tf_name, "_data.csv"))
# })


# Loop through each data frame in list_of_dfs
for(tf_name in names(list_of_dfs)) {
  current_df <- list_of_dfs[[tf_name]]
  
  # Check if the current data frame has more than 10 rows
  if(nrow(current_df) > 10) {
    
    # Rename columns as specified: "norm_avg" to "vivo", "beta_val" to "HEK293"
    # and correct the HEK293 values by dividing by 100
    current_df <- current_df %>%
      rename(vivo = norm_avg, HEK293 = beta_val) %>%
      mutate(HEK293 = HEK293 / 100)
    
    # Select only the relevant columns for the boxplot
    df_for_plot <- current_df %>%
      select(vivo, HEK293, c1_avg, c2_avg, c3_avg)
    
    # Convert data from wide to long format for ggplot2
    df_long <- pivot_longer(df_for_plot, 
                            cols = c(vivo, HEK293, c1_avg, c2_avg, c3_avg),
                            names_to = "Condition", 
                            values_to = "Value")
    
    # Calculate the number of CpG sites (rows in the dataframe)
    num_of_CpG_sites <- nrow(current_df)
    
    # Create the boxplot using ggplot2 with y-axis scale from 0 to 1
    p <- ggplot(df_long, aes(x = Condition, y = Value, fill = Condition)) +
      geom_boxplot() +
      labs(title = paste("Boxplot for", tf_name),
           subtitle = paste("Number of CpG sites investigated:", num_of_CpG_sites),
           x = "Condition", y = "Value") +
      scale_fill_discrete(name = "Condition") +
      theme_minimal() +
      scale_y_continuous(limits = c(0, 1)) # Set y-axis limits from 0 to 1
    
    # Display the plot
    print(p)
  }
}

