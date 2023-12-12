# Load necessary libraries
library(tidyr)

# Create a sample data frame
patient_data <- read.csv("DNAm_age314.csv")

patient_data <- patient_data[,c("File_Name", "Vital.Status", "Diagnosis.Age", "adj_DNAm_age")]
colnames(patient_data) <- c("File_Name", "Status", "age", "adj_DNAm_age")


# Convert patient status to numeric binary values
patient_data$status_numeric <- as.numeric(patient_data$Status == "Alive")

# Calculate the Pearson correlation coefficient between age and numeric status
cor_coefficient <- cor(patient_data$age, patient_data$status_numeric, method = "pearson")

# Print the correlation coefficient
cat("Pearson correlation coefficient between age and patient status: ", cor_coefficient)


# Perform the hypothesis test for the Pearson correlation coefficient
cor_test <- cor.test(patient_data$age, patient_data$status_numeric, method = "pearson")

# Print the results
cat("Test statistic (t): ", cor_test$statistic, "\n")
cat("Degrees of freedom (df): ", cor_test$parameter, "\n")
cat("p-value: ", cor_test$p.value, "\n")

# Determine if the correlation is significant at the 0.05 level
if (cor_test$p.value < 0.05) {
  cat("The correlation is significant at the 0.05 level.\n")
} else {
  cat("The correlation is not significant at the 0.05 level.\n")
}
