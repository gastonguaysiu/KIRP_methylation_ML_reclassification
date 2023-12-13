
**h6.R**

h6.R is used to visualize the heatmaps on the differentially methylated CpG sites in the samples and groupings.

Input file requirements:

 - sig_meth_DE.csv: Lists significantly differentially methylated CpG sites within gene promoter regions. Notably, this subset differs from those associated with the biomarkers for cluster organization.
 - filenames.csv: Used as input files for the analysis.
 - simp_heat2.csv: Average beta values for each significantly differentially methylated CpG site within promoter regions.
 - notes2.csv: Basic clinical information and grouping data for each sample.

 Output:

 - heatmap.pdf: Heatmap displaying ln(beta) for all CpG sites per sample, with annotations for sample grouping.
 - simp_heat.pdf: Heatmap of ln(beta) for CpG sites averaged according to clustering/grouping.
 - simp_heat2.pdf: Heatmap of ln(beta) for biomarker CpG sites averaged according to clustering/grouping.


**side_bar.R**

 The side_bar.R script was later remove in favour of bubble plots that could more aptly represent the information fo the GO analysis. This script was used to build the sidebar plots in KIRP/mval_analysis/output_conv1/

