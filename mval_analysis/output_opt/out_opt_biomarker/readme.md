**Machine Learning Analysis Results README**


**Overview**
This subfolder contains an analysis on the ML-determined biomarkers. The files found here give us a better idea of what the ML algorithm might have found to be important biomarkers in determining the survival outcome of patients

**Contents**

Folders:

- hyp3p: contains lists of of the hypermethylated and hypomethylated CpG associated genes for each cluster of patients


Files:

 - bavl_all2.csv: contains the sum of all the CpG probes that were determined as biomarkers (convA-E), along with their beta values for each sample
 - out.csv: is the subset of the 450k Illumna manifest files that only contains information of the probes that were isolated as opt. biomarkers
 - outB.csv: is a column reduction of the out.csv files, and contains a more concise set of information
 - small_sum.csv: Contains average methylated beta values, standard deviation, p-value, and p-adjusted value for the biomarker CpG sites in tumor clusters and normal tissue samples.
 - genomic_feat.csv contains information of the relative abundance in genomic location/features of the CpG sites found in the 450k Illumina manifest files, as well as the relative abundance of the biomarkers. i.e the location relitve to the nearest CpG island, and relative location to the nearest genes.
 - david.csv: the GO and KEGG results from a GO anlysis of the biomarker associated genes using the david platform.
 - small_genes.csv: a list of genes associated to the biomarkers according to the Illumnia 450k manifest files
 - gene_fun.csv: contains the biomarkers and associated genes, along with a little bit of information on whether the gene is associated with KIRP, RCC, or other cancers
 - fisher_stats.csv, fisher_val.csv: contains the information to see if the biomarkers would be considered a different set than the differentially methylated site based on their distribution in the genomic location/features


Scripts:

 - run_small.R: take the small_sum.csv as input and then determine the different hypermethylated and hypomethylated genes for the hyp3p folder
  - run_small2.R: take the small_sum.csv as input and then determine the different genes for each cluster
 - fisher_test.R: determines the fisher value for the fisher_stats.csv
