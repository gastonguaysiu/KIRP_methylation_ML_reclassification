# miRNA Classification Algorithm Overview

This repository hosts a series of R scripts and data files developed to classify a validation dataset based on centroid positions determined through an ML process involving miRNA data. which is based on on centroid positions determined through a series of convergences and an EM K-means cycle. Below is a detailed explanation of each script and its role in the process.

## Scripts and Their Roles

### Phases 1 to 3: Convergence Series

- **Script Name:** `alt_miRNA_f1.R`
- **Purpose:** This script is responsible for executing phases 1 through 3 of the algorithmic process. It systematically runs all 10 convergences in series, setting the foundation for the centroid determination in phase 4. Each convergence iterates through a predefined set of operations to refine the dataset's classification parameters progressively.

### Phase 4: Determining Centroid Positions

- **Script Name:** `alt_fixed_position.R`
- **Purpose:** In phase 4, this script takes over to finalize the centroid positions. These centroids are crucial as they are used to predict the classification of the validation dataset. The script employs a fixed position approach to ensure the centroids are accurately determined and consistent for predictive purposes.

### Hyperparameter Optimization and Dimensionality Reduction

- **Script Name:** `alt_opt_miRNA_f1.R`
- **Purpose:** Aimed at optimizing miRNA hyperparameters through an EM (Expectation-Maximization) K-means cycle to achieve a significant reduction in dimensionality. Despite the intentions, this method was unsuccessful for predictive purposes due to overfitting. This is like due to their miRNAs tendency to affect a multitude of genes inadvertently along with a relatively low number of hyperparameters within the context of genetic data.

## Challenges with miRNA Hyperparameters

The primary challenge in employing miRNA data for k-means clustering lies in the relatively low number of hyperparameters (1485) that can potentially be used to define the dimensions of the clustering, along with miRNAs tendency to affect a multitude of genes. When the hyperparameters of the convergences are summed across convergences 1-10, results in a highly dimensional space that is too general. A singular instance of convergence is too specific, which leads to overfitting.

## Running the Scripts

1. Ensure R is installed on your system.
2. Download or clone this repository.
3. Open and run the scripts in RStudio or any R environment in the following sequence:
   - `alt_miRNA_f1.R` for phases 1-3.
   - `alt_fixed_position.R` for phase 4.
   - (Optional) `alt_opt_miRNA_f1.R` to observe the attempted dimensionality reduction and its effects.


