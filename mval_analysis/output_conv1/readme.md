**Machine Learning Analysis Results README**


**Overview**
This repository contains the results from the first convergence of a Machine Learning (ML) algorithm, focusing on analyzing CpG-associated genes and ML-determined biomarkers. Although the results show potential, a better analysis can be made when comparing the initial (e0) and optimized (opt) convergences, highlighting the framework's potential; so, I won't discuss the results of this subfolder.

**Contents**

Folders:

- em_probes: Contains information on ML-determined biomarkers.
- cl1, cl2, cl3: Each folder holds lists of CpG-associated genes for respective clusters.

Files:

 - filenames.csv, kirp_2018_clinical_data.tsv: Used as input files for the analysis.
 - em_see.csv: An updating file representing the optimal subset of CpG sites for clustering.
 - na_best.csv, nb_best.csv: Current best scores for the subset of CpG sites used in k-means clustering.
 - clinical_cluster_data.csv: Averages clinical data for each cluster.
 - stat_summary.csv: Contains average methylated beta values, standard deviation, p-value, and p-adjusted value for each CpG site in tumour clusters and normal tissue samples.
 - sig_meth_DE.csv: Lists significantly differentially methylated CpG sites within gene promoter regions. Notably, this subset differs from those associated with the biomarkers for cluster organization.
 - notes2.csv: Basic clinical information and grouping data for each sample.
 - simp_heat2.csv: Average beta values for each significantly differentially methylated CpG site within promoter regions.

Visualizations:

 - elbow1.pdf: Elbow graph showing the ideal number of groups/clusters for biomarker CpG sites.
