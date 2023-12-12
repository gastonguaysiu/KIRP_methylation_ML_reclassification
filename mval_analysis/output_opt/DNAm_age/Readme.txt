Please note that I did not develop this script but rather used Steve Horvath's script that accompanied his "DNA methylation age of human tissues
and cell types" paper.

Horvath, S. (2013). DNA methylation age of human tissues and cell types. Genome Biology, 14(10), R115. https://doi.org/10.1186/gb-2013-14-10-r115

-----------------------------------------

To use the 450k Illumina data instead of the 27k data provided in the code, you'll need to make the following changes:

1. Replace the probe annotation file for 450k data:
Replace the file used to load probe annotation from 27k (13059_2013_3156_MOESM21_ESM.csv) to the file containing the 450k probe annotation. You would need to obtain or create this file.

R Copy code
probeAnnotation27k = read.csv("your_450k_probe_annotation_file.csv")

2. Update the probeAnnotation21kdatMethUsed:
This file (13059_2013_3156_MOESM22_ESM.csv) contains the 21k CpG sites used for the methylation age prediction. You may need to obtain or create a similar file with the corresponding CpG sites from the 450k Illumina data.

R Copy code
probeAnnotation21kdatMethUsed = read.csv("your_450k_corresponding_21k_CpG_sites_file.csv")

3. Read DNA methylation data:
Replace the input data file (13059_2013_3156_MOESM26_ESM.csv) with your 450k Illumina DNA methylation data file. Ensure the file is formatted similarly, with samples in columns and probes in rows.
 - My file named t26I.csv is missing due to the size of the file


R Copy code
dat0 = read.csv.sql("your_450k_methylation_data_file.csv")

Note that these changes assume the 450k data file is formatted similarly to the 27k data file, and the probe names are consistent between the annotation files and the methylation data file. Also, ensure that the 450k probe annotation file and the corresponding 21k CpG sites file contain the same columns and information as the original files.

