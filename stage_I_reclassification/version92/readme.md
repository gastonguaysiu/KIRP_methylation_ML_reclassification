

The **run.R** script will produce a variety of files needed for analysis and the development of figures. 
 - the stat_summary.csv contains the average methylated beta values for all the analyzed CpG sites for each k-means organized tumour cluster and the average beta value associated with the normal tissue samples. The table also contains the standard deviation, p-value, and p-adjusted value for each site
 - sig_meth_DE.csv contains all of the significantly differentially methylated CpG sites that are located within the gene promoter regions. Please note that this is not the same subset of CpG sites associated with the biomarkers used to organize the different clusters into the optimal groupings
 - notes2.csv contains the basic clinical information and grouping information associated with each sample
 - simp_heat2.csv contains the average beta values of each CpG site within the promoter region that was significantly differentially methylated.
 - cl*_dm_genes*.csv is the file that contains the differentially methylated genes associated with the 450k manifest files, which are NOT limited to the promoter regions

The **run2.R** script will produce clinical_cluster_data.csv: this file contains the clustering and clinical data for all samples. It therefore allows us to benchmark our categorization against clinical features.

The **mval_run_sum.R** script will produce mval_sum.csv: this file contains similar data as the stat_summary.csv file; however, it is derived from the M-values (instead of beta values). mval_sum.csv contains the average M-value of each cluster and control, the standard deviation, the p-value, and the p-adjusted values compared to the healthy tissue controls.

