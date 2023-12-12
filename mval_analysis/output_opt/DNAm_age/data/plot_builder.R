# Load necessary libraries
library(tidyverse)

# Create example data
data <- read.csv("DNAm_age314_status.csv")
colnames(data) <- c('group', 'age', 'adj_DNAm_age', 'adj_RE')


# Assuming your data is in a dataframe called 'data'
# The columns are 'group', 'age', 'adj_DNAm_age', 'adj_RE'

# Reshape data to long format
long_data <- data %>%
  gather(key = "Age_type", value = "Value", age, adj_DNAm_age) %>%
  unite(Group, group, Age_type, sep = "_")

# Calculate mean relative error (MRE) for each group
mre_data <- data %>%
  gather(key = "Age_type", value = "Value", age, adj_DNAm_age) %>%
  gather(key = "RE_type", value = "RE_value", adj_RE) %>%
  unite(Group, group, Age_type, sep = "_") %>%
  group_by(Group) %>%
  summarize(MRE = mean(RE_value))

# Create a whisker plot (box plot) using ggplot2
age_whisker_plot <- ggplot(long_data, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot() +
  geom_text(data = mre_data, aes(x = Group, y = -Inf, label = paste("MRE:", round(MRE, 2))), vjust = -1.5) +
  theme_minimal() +
  labs(title = "Age and Adj. DNAm Age Distribution among Groups",
       x = "Group",
       y = "Age Value",
       fill = "Group") +
  scale_fill_brewer(palette = "Pastel1")

# Display the plot
print(age_whisker_plot)

