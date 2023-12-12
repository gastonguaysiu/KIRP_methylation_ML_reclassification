library(tidyverse)

summary <- read.csv("outB.csv")
rownames(summary) <- summary[,1]

### create some tables to export to sqlite3, where we restrict to only significant probes
### i.e difference in beta expression of 33% and a p-adjusted of less than 0.01

###################   probe_to_gene.sh  ######################

manifest <- read_csv("headless.csv")

temp <- summary
temp <- temp %>%               
  separate_rows(UCSC_RefGene_Name, sep=";")
temp <- as.data.frame(unique(temp))
temp <- temp %>%               
  separate_rows(UCSC_RefGene_Group, sep=";")
temp <- as.data.frame(unique(temp))


alx <- manifest[,c("IlmnID", "Next_Base", "Genome_Build", "CHR", "MAPINFO", "Strand", "Probe_SNPs", "Probe_SNPs_10", "UCSC_RefGene_Name", "UCSC_RefGene_Group", 
                        "UCSC_CpG_Islands_Name", "Relation_to_UCSC_CpG_Island", "DMR", "Enhancer", "Regulatory_Feature_Name", "Regulatory_Feature_Group", "DHS")]
alx <- alx %>%               
  separate_rows(UCSC_RefGene_Name, sep=";")
alx <- as.data.frame(unique(alx))
alx <- alx %>%               
  separate_rows(UCSC_RefGene_Group, sep=";")
alx <- as.data.frame(unique(alx))

nx <- c(0, 0)
nx <- data.frame(nx)
rownames(nx) <- c("alx", "opt")
colnames(nx) <- c("probes")

nx["alx","probes"] <- nrow(manifest)
nx["opt","probes"] <- nrow(summary)

nx$body <- nrow(subset(alx, UCSC_RefGene_Group == "Body"))
nx[2,"body"] <- nrow(subset(temp, UCSC_RefGene_Group == "Body"))

nx$five_utr <- nrow(subset(alx, UCSC_RefGene_Group == "5'UTR"))
nx[2,"five_utr"] <- nrow(subset(temp, UCSC_RefGene_Group == "5'UTR"))
nx$three_utr <- nrow(subset(alx, UCSC_RefGene_Group == "3'UTR"))
nx[2,"three_utr"] <- nrow(subset(temp, UCSC_RefGene_Group == "3'UTR"))

nx$TSS200 <- nrow(subset(alx, UCSC_RefGene_Group == "TSS200"))
nx[2,"TSS200"] <- nrow(subset(temp, UCSC_RefGene_Group == "TSS200"))
nx$TSS1500 <- nrow(subset(alx, UCSC_RefGene_Group == "TSS1500"))
nx[2,"TSS1500"] <- nrow(subset(temp, UCSC_RefGene_Group == "TSS1500"))


nx$S_Shore <- nrow(subset(alx, Relation_to_UCSC_CpG_Island == "S_Shore"))
nx[2,"S_Shore"] <- nrow(subset(temp, Relation_to_UCSC_CpG_Island == "S_Shore"))
nx$S_Shelf <- nrow(subset(alx, Relation_to_UCSC_CpG_Island == "S_Shelf"))
nx[2,"S_Shelf"] <- nrow(subset(temp, Relation_to_UCSC_CpG_Island == "S_Shelf"))
nx$N_Shore <- nrow(subset(alx, Relation_to_UCSC_CpG_Island == "N_Shore"))
nx[2,"N_Shore"] <- nrow(subset(temp, Relation_to_UCSC_CpG_Island == "N_Shore"))
nx$N_Shelf <- nrow(subset(alx, Relation_to_UCSC_CpG_Island == "N_Shelf"))
nx[2,"N_Shelf"] <- nrow(subset(temp, Relation_to_UCSC_CpG_Island == "N_Shelf"))
nx$Island <- nrow(subset(alx, Relation_to_UCSC_CpG_Island == "Island"))
nx[2,"Island"] <- nrow(subset(temp, Relation_to_UCSC_CpG_Island == "Island"))

nx$CDMR <- nrow(subset(alx, DMR == "CDMR"))
nx[2,"CDMR"] <- nrow(subset(temp, DMR == "CDMR"))
nx$RDMR <- nrow(subset(alx, DMR == "RDMR"))
nx[2,"RDMR"] <- nrow(subset(temp, DMR == "RDMR"))

nx$enhancer <- nrow(subset(alx, Enhancer == "TRUE"))
nx[2,"enhancer"] <- nrow(subset(temp, Enhancer == "TRUE"))

nx$DHS<- nrow(subset(alx, DHS == "TRUE"))
nx[2,"DHS"] <- nrow(subset(temp, DHS == "TRUE"))


write.csv(nx,"genomic_feat.csv")
