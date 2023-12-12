
Additionally unique to each seed, representing the subset of probes for the optimized/final phase of the EM K-means algorithm are the following:

 - em_probes: this folder contains information for the ML determined biomarkers

 - cl1, cl2, cl3: these folders contain the lists for CpG associted genes for each of the


 - filenames.csv, kirp_2018_clinical_data.tsv: these files are used for input as analysis.
 - em_see.csv: is an updating file representing the optimum subset of the CpG site to use in the clustering.

 - na_best.csv, nb_best: represent the current best scores for the subset of the CpG site used in the k-means clustering.
 - clinical_cluster_data.csv: this averages out the clincal data for each cluster
 - stat_summary.csv: contains the average methylated beta values for all the analyzed CpG sites for each k-means organized tumour cluster and the average beta value associated with the normal tissue samples. The table also contains the standard deviation, p-value, and p-adjusted value for each site
 - sig_meth_DE.csv: contains all of the significantly differentially methylated CpG sites that are located within the gene promoter regions. Please note that this is not the same subset of CpG sites associated with the biomarkers used to organize the different clusters into the optimal groupings
 - notes2.csv: contains the basic clinical information and grouping information associated with each sample
 - simp_heat2.csv: contains the average beta values of each CpG site within the promoter region that was significantly differentially methylated.

heatmap.pdf: contains a heatmap for the ln(beta) for all the CpG sites for each sample, with annotations to indicate the grouping of which each sample belongs to
simp_heat.pdf: contains a heatmap for the ln(beta) for all the CpG sites that is averaged out accourding to the clustering/grouping
simp_heat2.pdf: contains a heatmap for the ln(beta) for the biomaker CpG sites that is averaged out accourding to the clustering/grouping
elbow1.pdf: elbow graph that displays the ideal number of groups/cluster for the biomarker CpG sites
km_survival.pdf: this is the KM survival plot for the clusters
Screenshot from 2023-02-01 14-07-08.png

sex.png: sex distribution for each cluster
fga.png: boxplot distribution of the fraction of genome altered groupped based on each cluster





