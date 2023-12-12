# Load ggplot2 and dplyr packages
library(ggplot2)
library(dplyr)

# Read the data and rename columns
df <- read.csv("start.csv")
colnames(df) <- c("Category", "Term", "Count", "Percent", "P_val", "Benjamini", "Convergence")

# Transform P_val to get -log(P_val) with a lower cap of 1e-10
df <- df %>%
  mutate(P_val = pmax(P_val, 1e-10),
         neg_log_P_val = -log10(P_val))

# Function to create plots
create_plot <- function(data_frame, title) {
  data_frame %>%
    arrange(desc(Convergence), desc(P_val)) %>%
    mutate(
      Term = factor(Term, levels = rev(unique(Term))),
      Convergence = factor(Convergence, levels = c("opt", "conv_E", "conv_D", "conv_C", "conv_B", "conv_A")) # Reversed order
    ) %>%
    ggplot(aes(x = Convergence, y = Term, size = neg_log_P_val, color = Category)) +
    geom_point(alpha = 0.6) +
    scale_color_manual(values = c("blue", "red", "green", "purple")) + # Adjust colors as needed
    labs(title = title,
         size = "-log P value",
         color = "Category",
         y = "Term",
         x = "Convergence Set") +
    theme_light() +
    theme(legend.position = "right", axis.text.x = element_text(angle = 90, hjust = 1))
}

# Plot excluding KEGG_PATHWAY
df_excluding_kegg <- filter(df, Category != "KEGG_PATHWAY")
plot_excluding_kegg <- create_plot(df_excluding_kegg, "Enrichment Analysis (Excluding KEGG_PATHWAY)")
print(plot_excluding_kegg)

# Plot including only KEGG_PATHWAY
df_only_kegg <- filter(df, Category == "KEGG_PATHWAY")
plot_only_kegg <- create_plot(df_only_kegg, "Enrichment Analysis (Only KEGG_PATHWAY)")
print(plot_only_kegg)
