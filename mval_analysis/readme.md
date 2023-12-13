**Heatmap & KM-survival plots**

We identified a CpG island methylator phenotype (CIMP) in the first cluster (cl1), marked by predominant hyper-methylation at CpG promoter sites. This cluster showed the poorest survival rates irrespective of whether or not the ML algorithm underwent estimation maximization to converge on a subset of CpG sites. However, after the convergence of a subset of CpG sites, the ML algorithm favoured more strongly grouping patients with expected poor survival outcomes together which can be seen in the following figures (refer to the scripts h6.R to build heatmaps and run.R/notes2.csv for the data on building the KM-survival plots)


clinical features

Changes were also observed in the distribution of clinical features in kidney renal papillary cell carcinoma (KIRP) tumours after optimizing the methylomic classification of patient samples using the EM k-means algorithm. The clinical data summary of the ML methylomic groupings was done using the run2.R script.

Gender Distribution Shift in Clusters: While many clinical features remained similar between the initial (e0) and optimized clusters, a notable shift was observed in the male-to-female patient ratio in each group. In the group with poorer survival outcomes, the proportion of female patients decreased significantly from 70% at e0 to 44.4% after optimization. This second value more accurately reflect female the male ratio that was observed in patients that succumb to cancer

Tumor Stage Distribution: Tumor staging, crucial for predicting cancer prognosis and survival, varies from T1 (localized cancer) to T3 (advanced or metastatic cancer). The distribution of tumour stages within each cluster changed after applying the EM algorithm. Notably, the group with the best survival outcome saw an increase in the ratio of patients in stage T1, rising from 68.18% to 74.07%. Conversely, in the worst outcome group, there was a shift away from stage T2, with the T2 ratio decreasing from 66.67% to 36.36%, T1 from increasing 0% to 18.18%, and T3 increasing from 33.33% to 45.45%. This shift more accurately represents the distribution seen in the varied tumour stages of the patients who  succumb to cancer.


Categorizing the patient broadly based on their methylation profiles indicated that differentially methylated CpG-associated genes il linked with different survival outcomes. As one might expect, the patients with the worst survival outcome had the largest number of differentially methylated CpG-associated genes linked with KIRP genes. After estimation maximation to refine CpG sites to an optimized subset to categorize patients, the differentially methylated profiles show a stronger correlation with the known KIRP genes (refer to script iso_genes.R creates the data table information).
