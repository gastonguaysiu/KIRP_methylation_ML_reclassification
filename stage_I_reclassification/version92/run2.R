library(tidyverse)

print("getting clinical data info for each cluster")

notes2 <- read.csv("notes2.csv")
notes2 <- notes2[,-1]
colnames(notes2) <- c("File_Name", "Case_ID", "Sample_ID", "lc_from_IPDD", "Days_to_death", "Vital_Status", "File_ID", "Sample_Type", "cluster")

cluster1 <- filter(notes2, cluster == 1, Sample_Type == "Primary Tumor")
cluster2 <- filter(notes2, cluster == 2, Sample_Type == "Primary Tumor")
cluster3 <- filter(notes2, cluster == 3, Sample_Type == "Primary Tumor")

normal <- filter(notes2, Sample_Type == "Solid Tissue Normal")

######  reducing the given to the assoc. sample ID  ######

clinical_data <- read.table(file = 'kirp_2018_clinical_data.tsv', sep = '\t', header = TRUE)

dead <- filter(notes2, Vital_Status == "Dead")
dead <- unique(as.data.frame(dead$Sample_ID))

normalc <- as.data.frame(normal$Sample_ID)
cl1c <- as.data.frame(cluster1$Sample_ID)
cl2c <- as.data.frame(cluster2$Sample_ID)
cl3c <- as.data.frame(cluster3$Sample_ID)

colnames(dead) <- c("sample_ID")
colnames(normalc) <- c("sample_ID")
colnames(cl1c) <- c("sample_ID")
colnames(cl2c) <- c("sample_ID")
colnames(cl3c) <- c("sample_ID")

dead <- as.data.frame(substr(dead$sample_ID, 1, nchar(dead$sample_ID)-1))
normalc <- as.data.frame(substr(normalc$sample_ID, 1, nchar(normalc$sample_ID)-1))
cl1c <- as.data.frame(substr(cl1c$sample_ID, 1, nchar(cl1c$sample_ID)-1))
cl2c <- as.data.frame(substr(cl2c$sample_ID, 1, nchar(cl2c$sample_ID)-1))
cl3c <- as.data.frame(substr(cl3c$sample_ID, 1, nchar(cl3c$sample_ID)-1))

colnames(dead) <- c("sample_ID")
colnames(normalc) <- c("sample_ID")
colnames(cl1c) <- c("sample_ID")
colnames(cl2c) <- c("sample_ID")
colnames(cl3c) <- c("sample_ID")

# create a minimal table to compare the results
# initializing the data frame of new copy results
clinical_sum <- c(0, 0, 0, 0, 0)
clinical_sum <- data.frame(clinical_sum)
rownames(clinical_sum) <- c("dead", "normal","cl1", "cl2", "cl3")
colnames(clinical_sum) <- c("avg_age")

cl_temp <- merge(x = dead, y = clinical_data, by.x="sample_ID", by.y="Sample.ID", all = FALSE)

######  fill out data for each row  ######
x = 0

while (x < 5) {
  x <- x + 1
  
  if (x == 1 ) {
    cl_temp <- merge(x = dead, y = clinical_data, by.x="sample_ID", by.y="Sample.ID", all = FALSE)
  } else if (x == 2 ) {
    cl_temp <- merge(x = normalc, y = clinical_data, by.x="sample_ID", by.y="Sample.ID", all = FALSE)
  } else if (x == 3 ) {
    cl_temp <- merge(x = cl1c, y = clinical_data, by.x="sample_ID", by.y="Sample.ID", all = FALSE)
  } else if (x == 4 ) {
    cl_temp <- merge(x = cl2c, y = clinical_data, by.x="sample_ID", by.y="Sample.ID", all = FALSE)
  } else {
    cl_temp <- merge(x = cl3c, y = clinical_data, by.x="sample_ID", by.y="Sample.ID", all = FALSE)
  }
  
  clinical_sum[x,"avg_age"] = mean(cl_temp$Diagnosis.Age)
  
  clinical_sum[x,"avg_Aneuploidy_Score"] <- mean(cl_temp$"Aneuploidy.Score", na.rm = TRUE)
  
  clinical_sum[x,"avg_Buffa.Hypoxia.Score"] <- mean(cl_temp$"Buffa.Hypoxia.Score", na.rm = TRUE)
  
  clinical_sum[x,"avg_Mutation.Count"] <- mean(cl_temp$"Mutation.Count", na.rm = TRUE)
  
  clinical_sum[x,"avg_Fraction.Genome.Altered"] <- mean(cl_temp$"Fraction.Genome.Altered", na.rm = TRUE)
  
  clinical_sum[x,"Cancer.Metastasis.Stage_mx"] <- length(which(cl_temp$"American.Joint.Committee.on.Cancer.Metastasis.Stage.Code" == "MX"))/nrow(cl_temp)
  clinical_sum[x,"Cancer.Metastasis.Stage_m0"] <- length(which(cl_temp$"American.Joint.Committee.on.Cancer.Metastasis.Stage.Code" == "M0"))/nrow(cl_temp)
  clinical_sum[x,"Cancer.Metastasis.Stage_m1"] <- length(which(cl_temp$"American.Joint.Committee.on.Cancer.Metastasis.Stage.Code" == "M1"))/nrow(cl_temp)
  
  clinical_sum[x,"Neoplasm.Disease.Lymph.Node.Stage_nx"] <- length(which(cl_temp$"Neoplasm.Disease.Lymph.Node.Stage.American.Joint.Committee.on.Cancer.Code" == "NX"))/nrow(cl_temp)
  clinical_sum[x,"Neoplasm.Disease.Lymph.Node.Stage_n0"] <- length(which(cl_temp$"Neoplasm.Disease.Lymph.Node.Stage.American.Joint.Committee.on.Cancer.Code" == "N0"))/nrow(cl_temp)
  clinical_sum[x,"Neoplasm.Disease.Lymph.Node.Stage_n1"] <- length(which(cl_temp$"Neoplasm.Disease.Lymph.Node.Stage.American.Joint.Committee.on.Cancer.Code" == "N1"))/nrow(cl_temp)
  clinical_sum[x,"Neoplasm.Disease.Lymph.Node.Stage_n2"] <- length(which(cl_temp$"Neoplasm.Disease.Lymph.Node.Stage.American.Joint.Committee.on.Cancer.Code" == "N2"))/nrow(cl_temp)
  
  clinical_sum[x,"Tumor.Stage_t1"] <- length(which(cl_temp$"American.Joint.Committee.on.Cancer.Tumor.Stage.Code" == "T1"))/nrow(cl_temp)
  clinical_sum[x,"Tumor.Stage_t1a"] <- length(which(cl_temp$"American.Joint.Committee.on.Cancer.Tumor.Stage.Code" == "T1A"))/nrow(cl_temp)
  clinical_sum[x,"Tumor.Stage_t1b"] <- length(which(cl_temp$"American.Joint.Committee.on.Cancer.Tumor.Stage.Code" == "T1B"))/nrow(cl_temp)
  clinical_sum[x,"Tumor.Stage_t2"] <- length(which(cl_temp$"American.Joint.Committee.on.Cancer.Tumor.Stage.Code" == "T2"))/nrow(cl_temp)
  clinical_sum[x,"Tumor.Stage_t2a"] <- length(which(cl_temp$"American.Joint.Committee.on.Cancer.Tumor.Stage.Code" == "T2A"))/nrow(cl_temp)
  clinical_sum[x,"Tumor.Stage_t2b"] <- length(which(cl_temp$"American.Joint.Committee.on.Cancer.Tumor.Stage.Code" == "T2B"))/nrow(cl_temp)
  clinical_sum[x,"Tumor.Stage_t3"] <- length(which(cl_temp$"American.Joint.Committee.on.Cancer.Tumor.Stage.Code" == "T3"))/nrow(cl_temp)
  clinical_sum[x,"Tumor.Stage_t3a"] <- length(which(cl_temp$"American.Joint.Committee.on.Cancer.Tumor.Stage.Code" == "T3A"))/nrow(cl_temp)
  clinical_sum[x,"Tumor.Stage_t3b"] <- length(which(cl_temp$"American.Joint.Committee.on.Cancer.Tumor.Stage.Code" == "T3B"))/nrow(cl_temp)
  
  # clinical_sum[x,"Primary.Lymph.Node.Presentation.Assessment_yes.over.no"] <- (length(which(cl_temp$"Primary.Lymph.Node.Presentation.Assessment" == "Yes")))/(length(which(cl_temp$"Primary.Lymph.Node.Presentation.Assessment" == "No")))
  clinical_sum[x,"Radiation.Therapy_yes.over.no"] <- (length(which(cl_temp$"Radiation.Therapy" == "Yes")))/(length(which(cl_temp$"Radiation.Therapy" == "No")))
  clinical_sum[x,"avg_Ragnum.Hypoxia.Score"] <- mean(cl_temp$"Ragnum.Hypoxia.Score", na.rm = TRUE)
  clinical_sum[x,"Sex_male_over_female"] <- (length(which(cl_temp$"Sex" == "Male")))/(length(which(cl_temp$"Sex" == "Female")))
  
}

write.csv(clinical_sum,"clinical_cluster_data.csv")