
**mval_run_sum.R**

This script is used to build the small_sum.csv, which focuses on the biomarkers used for grouping the patient samples.

Input file requirements:

 - beta_val_all.csv: contains all the CpG Beta probe values for all the samples
 - mval_all.csv: contains all the CpG Beta probe values for all the samples
 - em_see.csv: contain a list of all the CpG site, along with an index of which should be used as biomarkers in clustering for the k-means algorithm
 - filenames.csv: Used as input files for the analysis.
 - ref4.csv: basic clinical information of each sample
 - headless.csv: contains the 450k Illumina manifest table, which is important information about the CpG probes

**run_sum.R**

Input file requirements: Same as the mval_run_sum.R

 Output:

 - small_sum.csv: Contains average methylated beta values, standard deviation, p-value, and p-adjusted value for the biomarker CpG sites in tumour clusters and normal tissue samples.
  - notes2.csv: Basic clinical information and grouping data for each sample.
 - sig_probes.csv: contains the list of CpG probes that were used as biomarkers for the ML algorithm for the grouping of the samples
 - conv_genes.csv: the genes that are associated with the CpG biomarkers according to the 450k Illumina manifest file

