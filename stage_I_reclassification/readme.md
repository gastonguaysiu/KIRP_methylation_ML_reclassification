# ML-reclassification
ML-KIRP reclassification using k-means and EM

Estimation Maximization (EM) is a general-purpose iterative algorithm used in machine learning to maximize a posteriori (MAP) estimates of parameters in statistical models; in our case, the parameter is the probes used in the k-means clustering. Given that we know the outcome of the patients, the EM part of the algorithm would be considered supervised, and the K-means would be regarded as unsupervised. EM is beneficial in situations where some of the variables or parameters are unobserved or hidden. The idea behind EM is to alternate between estimating the missing or hidden variables and updating the model parameters based on these estimates. This process continues until the estimates of the model parameters converge to a maximum or a local maximum.


**Folder table of contents**

- the version73_conv_alt_seed folder contains four more subfolders, each containing a different trail. The clust73em_bval subfolder contains a trial attempted on the beta values rather than M-values for each sample. The other three sub-folders contain the same trials attempted using a different random number seed to evaluate the reproducibility of the data.
- The version74_opt folder contains data for phase 3 of the  clust73em_bval and clust73em_mval trails. I.E it contains the CpG subset results of the EM k-means algorithm
- the version92 folder is the first step in processing the results from the methylation M-values for further analysis.

MISSING FILES: The scripts will not run as they are due to the file size restrictions. Users who wish to use the scripts should download methylation data from https://gdc.cancer.gov/ or other repositories. And convert them to M values. Additionally, the 450k Illumina manifest file in csv format is required for some of the scripts to run. 

**The differential methylation EM K-means workflow comprises the following steps.**

Phases 1-2 of the EM k-means algorithm use the **pip4.R** script. Each new iteration requires a new em_see.csv, mval_all.csv and bval_all.csv files, the unused/non-selected subset of CpG probes from the previous iteration. Running the **alt_see.R** script builds these new files for each iteration.
Phase 3 requires the user to manually create a new em_see.csv file based on the previous iterations (Conv. A - Conv. E).

Phase 1, called Initialization (e0), involves using all the available probes as a baseline. This means that the (probe ✕ sample) matrix initially contains a row for each DNA CpG probe with no NA value in any samples. The matrix is then grouped into clusters using k-means clustering. Each cluster is scored based on clinical data and organized into predicted survival outcome groups.

Phase 2, called the Estimation Maximization Inner Loop, involves copying the probes from the current best list and creating a new estimation by randomly adding and removing a few probes. 
1. The probes from the current best list are copied to create a trial clone.
2. The trial clone randomly adds and removes a unique set of probes to create a new estimation. Up to half of the current probes are removed, and up to the same number of unused probes are added.
3. Unsupervised K-means clustering is performed on the samples using the probes in the new estimation, resulting in three groups (k = 3).
4. Scoring uses clinical data to compare the new estimations to previous ones. A new best-estimated probe list is saved if one of the following criteria is met (in order of priority):
 - The clustering of the worst outcome group is more accurate.
 - The best OS group has a better censored-to-death ratio without hindering the worse outcome group.
 - Scoring improves because more samples are added to the extreme ends of the groups, indicating that patients are correctly classified into the best and worst predicted survival groups without losing accuracy in prediction.
 - The new estimation has fewer probes, but the scoring remains the same. This maximizes the weight of new probes or removes probes in the k-means clustering step.
5. The best estimation is fed back into the k-means clustering in step A. The process is repeated until a predetermined number of new estimation trials fail to improve the best estimation, indicating that the best list of probes for k-means clustering to segregate survival outcomes in KIRP tumour samples has been reached.


Phase 3, called Estimation Maximization Optimization, involves saving the list of probes that converged on the best results in Phase 2.5 and removing them from the pool of possible probes to use. The algorithm loops back to Phases 1 and 2, slightly reducing the number of probes used. This process is repeated for five iterations (Conv.A - Conv.E). Finally, the list of probes from the instances is initialized into the e0 step before being fed back into the algorithm in Phase 1, resulting in an optimized list of probes.


![flowchar](https://github.com/gastonguaysiu/KIRP/blob/main/stage_I_reclassification/ML_flow_KIRP.png?raw=true)



**The formulation of the scoring in the EM algorithm ==> definition of fitness**

Patients with cancer either die or discontinue checkups, but stopping checkups doesn't always imply survival. The time since diagnosis (lc_from_IPDD) and the date of death post-diagnosis are crucial data points. The F, calculated as

F=D−(l/α)

Where 'D' is the death count, 'l' is days since lc_from_IPDD, and 'α' is the median survival time, identifying clusters with the lowest survival rates. Conversely, the R, given by

R=(D/L)×(d/β),

where 'D' is the death count, 'L' is the sample size, 'd' survival days, and 'β' is the median lc_from_IPDD, predicts the best survival outcomes. Zero is avoided in calculations by assigning 0.01 to 'D' and 'd' when no deaths are recorded.
The algorithm's efficiency decreases with more data, so we focus on fewer significant CpG probes. A set of probes is optimal if smaller yet achieves the same score, balancing computational efficiency with accuracy.

