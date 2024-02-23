# EM Subfolder Overview

This `em` subfolder contains data and scripts associated with the miRNA classification algorithm, focusing on the Expectation-Maximization (EM) K-means cycle utilized during the algorithmic process. Within this subfolder, you will find two key subfolders: `non_norm_conv_1-10` and `norm_conv0-10`, each containing data processed under different conditions for the miRNA data.

## Subfolders and Their Data

### `non_norm_conv_1-10` Subfolder

The `non_norm_conv_1-10` subfolder houses data files associated with non-normalized miRNA data. These data files were the product of phases 2 and 3 of the machine learning cycle. The raw values used are based on reads per million miRNA mapped (RPM). This approach maintains the original distribution of the data, providing a unique perspective on the miRNA profiles without the influence of normalization procedures.

### `norm_conv0-10` Subfolder

Conversely, the `norm_conv0-10` subfolder contains data that has been normalized using the reads RPM methodology, applied row-wise. This normalization ensures that the data is scaled appropriately across different samples, allowing for a more standardized comparison of miRNA expression levels.

## Importance of Different Data Processing Approaches

The inclusion of both normalized and non-normalized data sets allows for a comprehensive analysis of the miRNA classification process. By comparing the outcomes of each approach, researchers can better understand the impact of data normalization on the clustering and classification results of the miRNA data, providing insights into the most effective methods for miRNA-based data interpretation.

## Utilizing the Data

Researchers and analysts can delve into these subfolders to explore the variations in miRNA classification outcomes based on the data processing techniques applied. It is crucial to understand the context and methodology behind each dataset to accurately interpret the results and leverage them for further research or validation purposes.

Remember to follow the specific script execution order mentioned in the main `readme.md` to ensure the proper application of these data sets within the miRNA classification algorithm framework.
