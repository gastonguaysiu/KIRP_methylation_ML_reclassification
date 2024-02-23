# Hybrid Predictive Model Implementation

## Overview
This document details implementing a hybrid predictive model that employs machine learning (ML) techniques to analyze differential methylation data. The model's core revolves around the strategic analysis of CpG probe sites for predictive accuracy in patient survival outcomes by implementing a combination of estimation-maximization and k-means clustering.


## Normalization of miRNA and RNAseq data
A log transformation is applied to the TPM or RPM values to mitigate the influence of highly expressed genes, adding a small constant (0.0001) to avoid undefined values for zeros. This step makes the data more normally distributed, a common requirement for many statistical analyses. Then, the data is centred and scaled (z-score normalization) using the formula (x - mean(x)) / sd(x) applied row-wise. This step ensures that each gene's expression is measured in terms of its deviation from its mean expression level, making comparisons across genes more meaningful. Normalization is performed on a per-transcript/miRNA basis across all samples.

## Data Handling
The dataset is partitioned into two segments:
- **Training Data**: Constitutes 90% of the samples, selected randomly.
- **Testing Data**: Comprises the remaining 10%, which is set aside for model validation.

## Predictive ML Workflow
The workflow is segmented into phases, each designed to refine the model's accuracy through a series of iterations. A detailed flowchart (Figure 2) illustrates this process.

![flowchar](https://github.com/gastonguaysiu/KIRP/blob/main/stage_II_prediction/ML_lock.png?raw=true)

### Phase 1: Initialization (e0)
- Utilizes all available probes to form a comprehensive (probe ✕ sample) matrix, ensuring no missing values across samples.
- K-means clustering groups the matrix into distinct populations/clusters based on differential methylation patterns, scored against clinical outcomes to predict survival rates.

### Phase 2: Estimation Maximization Inner Loop
- Proceeds with the best probe list from Phase 1, generating trial clones by randomly adding or removing probes to refine the estimation.
- This dynamic adjustment involves unsupervised K-means clustering, targeting three main groups (k = 3), with the aim of enhancing survival outcome prediction.
- Clinical data scoring compares new estimations against previous ones, focusing on accuracy in worst outcome group clustering, optimizing the censored-to-death ratio, and ensuring precise classification into survival groups.
- The process iterates until a peak in estimation quality is identified, marking the optimal list of probes for segregating survival outcomes in the studied samples.

### Phase 3: Enhanced Model Robustness
- Similar to the Estimation Maximization Inner Loop, but emphasizes model robustness over mere optimization.
- Prioritizes reintegrating probes when scoring matches the best estimation, enhancing the model's ability to predict with higher accuracy and robustness.
- This phase also includes the Estimation Maximization Optimization step, which conserves the probe list yielding the best results and iterates the process with a reduced probe set across ten iterations (Conv.1 - Conv.10), culminating in an optimized probe list.

### Phase 4: Final Optimization
- Compiles the refined list of probes from previous iterations to determine centroid locations in hyperdimensional space using K-means clustering.
- The centroid locations are saved, providing a stable basis for future sample group predictions based on CpG probe analysis.

## Scoring Metric
The F1 score to categorize patients into groups with the worst overall survival outcomes, is used for evaluating identifying CpG sites used as hyperparameters. The F1 score is a balance between precision and recall, crucial for handling imbalanced data classes in binary classification problems.

### Calculations
- **Precision** " = D / (D + (l/α)) "
- **Recall** " = D / (D + σ) "
- **F1 Score** " = 2 * (precision * recall) / (precision + recall) "

Where:
- D = Count of patients within the cluster who died due to cancer.
- l = Number of days counted in lc_from_IPDD.
- α = Median survival period before patients died of cancer across all patient data.
- σ = Number of known patients that succumbed to cancer in another cluster.

## Conclusion
This hybrid model integrates k-means clustering with estimation maximization, aiming for robust optimization in predicting patient groupings based on CpG probe analysis. The use of F1 scores for hyperparameter evaluation emphasizes the model's effectiveness in binary classification problems, particularly in scenarios with imbalanced datasets.
