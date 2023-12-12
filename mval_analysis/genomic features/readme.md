**Fisher's Exact Test Analysis Project**

**Introduction**

This repository contains R scripts for performing Fisher's Exact Test on two datasets: CpG Islands and Gene Groups. The analysis aims to explore the association between our optimized subset of biomarkers and the 450k manifest file using Fisher's Exact Test. This test is especially suitable for small sample sizes or when the assumptions of the Chi-squared test are not met.

Requirements

- R programming language
- R packages:
- tidyverse
- stats

Installation: Ensure you have R installed on your system. You can download it from CRAN. Once R is installed, open an R console and install the required packages by running:

**Project Structure**

CpG_Island.csv: Data file for CpG probes locations relative to the CpG Island analysis.

RefGene_Group.csv: Data file for CpG probes locations relative to the nearest genes analysis.

fishers_test_script.R: R script for conducting Fisher's Exact Test on the datasets.

**Running the Analysis**
To run the analysis, follow these steps:

- Clone this repository to your local machine.
- Open the fishers_test_script.R script in an R environment (like RStudio).
- Set the working directory to the location of the cloned repository.
- Run the script in the R environment.

The script performs Fisher's Exact Test on CpG Island and Gene Group datasets and calculates the p-values using simulated p-values for large sample sizes.

**Results**

The Fisher's Exact Test results indicate p-values that are statistically significant for the CpG_Island.csv for the comparison of the biomarkers to the 450k manifest file. The results are stored in the following variables:

CpG_island_p: 0.04973

CpG_gene_p: 0.1234
