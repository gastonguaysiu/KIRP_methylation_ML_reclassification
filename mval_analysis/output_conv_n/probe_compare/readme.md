Folders:

- fig_convN: contains dot plots and volcano plots that combine all the isolated biomarkers of the different convergences of the different runs of the ML algorithm.
-  fig_opt contains the dot plot and volcano plots that compare the change of the M-values to the controls scaled according to the standard deviation. These figures demonstrate that the biomarkers have a larger change in M-value than would be expected in a random sampling of CpG sites and that as the expected outcome of the cluster worsens, the M-value change will increase

Description of files:

 - mval_opt.csv: This file contains details of CpG sites' biomarkers used in the optimized clustering of the samples. Rows are the probes; columns include average M value, standard deviation, p-value, and adjusted p-value.  
 - mval_opt2.csv: this is just a reduced list of CpG site biomarker details. This file excludes CpG sites with p-values that are too high or don't have enough change between the tumour grouping and control. 
 - mval_optB.csv: a table that contains fewer rows; instead of having a row for each grouping/cluster, the 


 - mval_opt_sum.csv: this is the table that contains the average M value, standard deviation, p-value, and adjusted p-value for each cluster on all the CpG sites analyzed, run using the optimized biomarkers for k-means clustering    
 - mval_sum.csv: this contains the same information as mval_opt.csv but includes the information of the biomarkers used when a different random seed for the ML algorithm was used --- this file is not used to build any figure

- bval_opt2.csv, bval_opt.csv, bval_sum2.csv, bval_sum.csv: contains the same information as the mval* file but uses the beta values instead.

- used_opt.csv: This is a list of all the CpG sites that could be used for k-means clustering in the ML algorithm; the e0 = 1 rows are the biomarkers used in the optimized ML algorithm

- dot_plot.R: some preliminary script to visualize the biomarkers in reference to the other CpG sites as either dot plots or volcano plots, meant to answer the question, are these sites different by p-value, or change compared to the control or the other clusters?
