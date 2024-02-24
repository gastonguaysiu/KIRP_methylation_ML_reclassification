# RNA_seq Classification Algorithm Overview

This repository hosts a series of R scripts and data files developed to classify a validation dataset based on centroid positions determined through an ML process involving RNA_seq data, which is based on centroid positions determined through a series of convergences and an EM K-means cycle. Below is a detailed explanation of each script's role in the process.

## Scripts and Their Roles within the EM subfolder

### Hyperparameter Optimization and Dimensionality Reduction

- **Script Name:** `split_data.R`
- **Purpose:** Aimed at splitting the initial data set into the trainning (90%) and validation (10%) dataset. One key note about this is that the name of the files need to be cross reference across the different data types (CpG450K, or miRNA), should one choose to combine the data for a multi-omic predictive evaluation.

### Phases 1 to 3: Convergence Series

- **Script Name:** `alt_RNAseq_f1.R`
- **Purpose:** This script executes phases 1 through 3 of the algorithmic process. It systematically runs all 10 convergences in series, setting the foundation for the centroid determination in phase 4. Each convergence iterates through a predefined set of operations to progressively refine the dataset's classification parameters.

### Phase 4: Determining Centroid Positions

- **Script Name:** `alt_fixed_position.R`
- **Purpose:** In phase 4, this script takes over to finalize the centroid positions. These centroids are crucial as they are used to predict the classification of the validation dataset. The script employs a fixed position approach to ensure the centroids are accurately determined and consistent for predictive purposes.

### Hyperparameter Optimization and Dimensionality Reduction

- **Script Name:** `alt_opt_RNAseq_f1.R`
- **Purpose:** Aimed at optimizing RNA_seq hyperparameters through an EM (Expectation-Maximization) K-means cycle to achieve a significant dimensionality reduction. Despite the intentions, this method was unsuccessful for predictive purposes due to overfitting. This is due to their RNA_seq's tendency to affect a multitude of genes inadvertently, along with a relatively low number of hyperparameters within the context of genetic data.

## Challenges with RNA_seq Hyperparameters

The primary challenge in employing RNA_seq data for k-means clustering lies in the relatively low quality of the data brought on by the purity of the tumour. In the RNA_seq data, the normal tissue sample chare the highest similarity to the worst OS group based on PCA analysis.

## Running the Scripts

1. Ensure R is installed on your system.
2. Download or clone this repository.
3. Open and run the scripts in RStudio or any R environment in the following sequence:
   - `alt_RNAseq_f1.R` for phases 1-3.
   - `alt_fixed_position.R` for phase 4.
   - (Optional) `alt_opt_RNAseq_f1.R` to observe the attempted dimensionality reduction and its effects.

## MISSING FILES from RNA_seq data

trainingData.csv --- 274.9MB  --- 8 June 2023
The file trainingData.csv is a subset database for the normalized TPM values associated with all transcripts (including different isoforms) from patients that were use to train the ML model

tpm_RNA_seq2.csv --- 93.2MB  --- 10 June 2023
This file is a comprehensive database for the raw TPM transcript values associated collected from all patients. Contains the TPM value for all the patients.

tpm_mod_RNA_seq.csv --- 94.8MB  --- 6 June 2023
This file is a comprehensive database for the raw TPM transcript values associated collected from all patients; along with tow addition columns, on for the gene name associated with the transcript, and on for the type of transcript (e.g. protein coding)
