# RNA_seq Classification Algorithm Overview

This repository hosts a series of R scripts and data files developed to classify a validation dataset based on centroid positions determined through an ML process involving RNA_seq data, which is based on centroid positions determined through a series of convergences and an EM K-means cycle. Below is a detailed explanation of each script's role in the process.

## Scripts and Their Roles within the EM subfolder

### PCA for Dimensionality Reduction

- **Script Name:** `pca2.R`, `pca3.R`
- **Purpose:** The PCA analysis was to determine the similarity of the 3 clustered to the normal samples. The worst OS group was similar to the normal tissue samples, which was contrary to the CPG data. This is likely due to low tumour purity with the RNA_seq data.


## Running the Scripts (optional)

1. Ensure R is installed on your system.
2. Download or clone this repository.
3. Build out the missing files
4. Open and run the scripts in RStudio or any R environment in the following sequence:
   - `pca2.R`

## MISSING FILES from RNA_seq data

tpm_RNA_seq2.csv --- 93.2MB  --- 10 June 2023
This file is a comprehensive database for the raw TPM transcript values associated collected from all patients. Contains the TPM value for all the patients.

tpm_mod_RNA_seq.csv --- 94.8MB  --- 6 June 2023
This file is a comprehensive database for the raw TPM transcript values associated collected from all patients; along with tow addition columns, on for the gene name associated with the transcript, and on for the type of transcript (e.g. protein coding)
