**Machine Learning Analysis Results README**


**Overview**
This repository contains the results from the first convergence of a Machine Learning (ML) algorithm, focusing on analyzing CpG-associated genes and ML-determined biomarkers. Analysis on the efficacy of using my framework in the methylomic clasification of disease can be made when comparing the initial (e0) and optimized (opt) convergences. As such many of the comparisons between methylation profiles, and their implications are found within this subfolder.

**Contents**

Folders:

- diff_meth_genes: contains CpG-associated genes lists for respective clusters.
- DNAm_age: contains information and figure for an analysis done on the age of the optimized methylomicly classified patient compare to their clinical counterpart
- enrichment: contains the the genes and associated GO of the methylomicaly calsiffied groups, and the overlap found in the GO analysis for each grouping
- gene_counts: is a continuation of the enrichment analysis, with brings in the lsit of KIRP associated gene found by ncbi
- out_opt_biomarker is an extended analysis on just the biomarkers used in optimizing the methylomic classification of the patients
- pdfs: contains some of the Visualizations of the data, such as heatmaps, and elbow graph


Files:

 - filenames.csv, kirp_2018_clinical_data.tsv: Used as input files for the analysis.
 - em_see.csv: An updating file representing the optimal subset of CpG sites for clustering.
 - na_best.csv, nb_best.csv: Current best scores for the subset of CpG sites used in k-means clustering.
 - clinical_cluster_data.csv: Averages clinical data for each cluster.
 - stat_summary.csv: Contains average methylated beta values, standard deviation, p-value, and p-adjusted value for each CpG site in tumor clusters and normal tissue samples.
 - bval_summary_promo.csv: Contains average methylated beta values, standard deviation, p-value, and p-adjusted value for all the promoter associated CpG site in tumor clusters and normal tissue samples.
 - notes.csv, notes2.csv: Basic clinical information and grouping data for each sample.
 - simp_heat2.csv: Average beta values for each significantly differentially methylated CpG site within promoter regions.
 - clinical_cluster_data_mval_opt.ods: contains the clinical data for the e0 and opt methylomic k-means determined grouping for comparison.

Visualizations:

 - km_survival.png: Kaplan-Meier survival plot for the clusters.
 - Screenshot from 2023-02-04 13-22-20.png
 - sex.png: Sex distribution across each cluster.
 - fga.png: Boxplot showing the distribution of the altered genome fraction based on each cluster.
