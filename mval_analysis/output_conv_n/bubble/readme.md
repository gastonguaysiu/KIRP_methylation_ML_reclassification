Each of the enrichment*.csv files contains a table of the GO and KEGG information obtained by the Database for Annotation, Visualization and Integrated Discovery (DAVID) using the lists of genes associated with the CpG biomarkers.

The start.csv file contains the combined information of all the other tables and an additional column to indicate the source.

**bubble2.R**: builds a bubble plot for the GO and KEGG enrichment analysis using the start.csv file. The script is designed to run in Rstudio, and the visualizations must be exported manually.
