**Machine Learning Analysis Results README**


**Overview**
This repository contains the results from the first convergence of a Machine Learning (ML) algorithm, focusing on analyzing CpG-associated genes and ML-determined biomarkers. Analysis on the efficacy of using my framework in the methylomic clasification of disease can be made when comparing the initial (e0) and optimized (opt) convergences.

**Contents**

Folders:

- gene_count: contains CpG-associated genes lists for respective clusters.

Files:

 - kirp_2018_clinical_data.tsv: Used as input files for the analysis.
 - em_see.csv: An updating file representing the optimal subset of CpG sites for clustering.
 - na_best.csv, nb_best.csv: Current best scores for the subset of CpG sites used in k-means clustering.
 - clinical_cluster_data.csv: Averages clinical data for each cluster.
 - stat_summary.csv: Contains average methylated beta values, standard deviation, p-value, and p-adjusted value for each CpG site in tumor clusters and normal tissue samples.
 - notes2.csv: Basic clinical information and grouping data for each sample.
 - simp_heat2.csv: Average beta values for each significantly differentially methylated CpG site within promoter regions.

Visualizations:

 - heatmap.pdf: Heatmap displaying ln(beta) for all CpG sites per sample, with annotations for sample grouping.
 - simp_heat.pdf: Heatmap of ln(beta) for CpG sites averaged according to clustering/grouping.
 - simp_heat2.pdf: Heatmap of ln(beta) for biomarker CpG sites averaged according to clustering/grouping.
 - elbow1.pdf: Elbow graph showing the ideal number of groups/clusters for biomarker CpG sites.
 - km_survival.png: Kaplan-Meier survival plot for the clusters.
 - Screenshot from 2023-02-01 14-05-37.png
 - sex.png: Sex distribution across each cluster.
 - fga.png: Boxplot showing the distribution of the altered genome fraction based on each cluster.
