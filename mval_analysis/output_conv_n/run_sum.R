library(tidyverse)

# em_see file contains a matrix of i rows each to represent a probe.
# em_see contains 2 rows, a mater list of all probes, and the e0/best set of probes
# first column contains all the probe ID
# the second column contains 1s on probe ID for the current best estimation
em_see <- read.csv("em_see.csv")
em_see <- em_see[,-1]

em_see[em_see == 0] <- NA
new <- na.omit(em_see)
new <- as.data.frame(new[,1])
colnames(new) <- "Name"

# The beta and mval files contain a table each;
# beta-val contains a table of i probes and j samples, each cell representing the beta score
# where as mval contain cells representing the logit transformed beta value
# the merge function acts like a inner join in this case
bval_all <- read_csv("beta_val_all.csv")
mval_all <- read_csv("mval_all.csv")

beta_ex <- merge(x = new, y = bval_all, by.x = "Name", by.y = "...1")
mx <- merge(x = new, y = mval_all, by.x = "Name", by.y = "...1")

# pre-processing the files to remove any unwanted column taken up on 
# the sql inner join and re-enter corrected files-names for column names
# then transposing the matrix
rownames(mx) = mx[,1]
mx <- mx[,-1]
filenames <- read.csv("filenames.csv")
colnames(mx) <- t(filenames)

set.seed(2)
clusters <- as.data.frame(filenames)
clusters$cluster <- kmeans(t(mx),3,iter.max=10000)$cluster

# create a table with important information about
# each sample, such as file name, case id, sample id, vital statue, 
# and associate cluster
ref <- read.csv("ref4.csv")
notes <- merge(x = ref, y = clusters, by.x = "File_Name", by.y = "names")
notes1 <- notes[,-c(4, 5, 9, 10, 11, 13, 14, 15)]

#create some table with the info for the cluster tumor samples and another set 
# of tables with cluster normal samples.
cluster1 <- filter(notes1, cluster == 1, Sample_Type == "Primary Tumor")
cluster2 <- filter(notes1, cluster == 2, Sample_Type == "Primary Tumor")
cluster3 <- filter(notes1, cluster == 3, Sample_Type == "Primary Tumor")

normal <- filter(notes1, Sample_Type == "Solid Tissue Normal")
normal1 <- filter(notes1, cluster == 1, Sample_Type == "Solid Tissue Normal")
normal2 <- filter(notes1, cluster == 2, Sample_Type == "Solid Tissue Normal")
normal3 <- filter(notes1, cluster == 3, Sample_Type == "Solid Tissue Normal")

###################   stats.R   ######################

print ("pre-processing and loading up beta_current...")

beta_current <- beta_ex

temp <- beta_current[,1]
temp <- as.list(t(temp))
beta_current <- beta_current[,-c(1,2)]
rownames(beta_current) <- c(temp)
colnames(beta_current) <- t(filenames)

beta_t <- t(beta_current)
temp <- as.data.frame(rownames(beta_t))
colnames(temp) <- "file_name"
beta_t <- cbind(temp, beta_t)

### seperating information into data frames to later get statistical info (pre-processing) ###
print ("building a data frame for normal + each cluster")

normalB <- as.data.frame(normal$File_Name)
colnames(normalB) <- c("File_Name")
normalB <- merge(x=normalB, y=beta_t, by.x="File_Name", by.y="file_name", all = FALSE)
temp <- as.list(t(normalB[,1]))
rownames(normalB) <- c(temp)
normalB <- t(normalB[,-1])

cluster1B <- as.data.frame(cluster1$File_Name)
colnames(cluster1B) <- c("File_Name")
cluster1B <- merge(x=cluster1B, y=beta_t, by.x="File_Name", by.y="file_name", all = FALSE)
temp <- as.list(t(cluster1B[,1]))
rownames(cluster1B) <- c(temp)
cluster1B <- t(cluster1B[,-1])

cluster2B <- as.data.frame(cluster2$File_Name)
colnames(cluster2B) <- c("File_Name")
cluster2B <- merge(x=cluster2B, y=beta_t, by.x="File_Name", by.y="file_name", all = FALSE)
temp <- as.list(t(cluster2B[,1]))
rownames(cluster2B) <- c(temp)
cluster2B <- t(cluster2B[,-1])

cluster3B <- as.data.frame(cluster3$File_Name)
colnames(cluster3B) <- c("File_Name")
cluster3B <- merge(x=cluster3B, y=beta_t, by.x="File_Name", by.y="file_name", all = FALSE)
temp <- as.list(t(cluster3B[,1]))
rownames(cluster3B) <- c(temp)
cluster3B <- t(cluster3B[,-1])

# create a function to be used in finding the standard deviation... fuck r base ##
rowVars <- function(x, na.rm=F) {
  # Vectorised version of variance filter
  rowSums((x - rowMeans(x, na.rm=na.rm))^2, na.rm=na.rm) / (ncol(x) - 1)
}

### finally getting those stats ###
print ("getting stats for norm + clusters...")

summary <- as.data.frame(rowMeans(normalB))
colnames(summary) <- ("norm_avg")
summary$norm_std <- sqrt(rowVars(normalB))
summary$c1_avg <- rowMeans(cluster1B)
summary$c1_std <- sqrt(rowVars(cluster1B))
summary$c2_avg <- rowMeans(cluster2B)
summary$c2_std <- sqrt(rowVars(cluster2B))
summary$c3_avg <- rowMeans(cluster3B)
summary$c3_std <- sqrt(rowVars(cluster3B))

### use Kolmogorov–Smirnov test of two samples (cluster vs norm) to get a p-value. KS tests shape of distribution unlike the t.test counterpart

summary$pval_norm_vs_c1 <- sapply(1:nrow(normalB), function(i) ks.test(as.vector(normalB[i,]), as.vector(cluster1B[i,]))$p)
summary$pval_norm_vs_c2 <- sapply(1:nrow(normalB), function(i) ks.test(as.vector(normalB[i,]), as.vector(cluster2B[i,]))$p)
summary$pval_norm_vs_c3 <- sapply(1:nrow(normalB), function(i) ks.test(as.vector(normalB[i,]), as.vector(cluster3B[i,]))$p)


# use Benjamini & Hochberg (1995) fdr, because stat quest youtube video explained this particular fdr

summary$padj_norm_vs_c1 <- p.adjust(summary$pval_norm_vs_c1, method = "BH")
summary$padj_norm_vs_c2 <- p.adjust(summary$pval_norm_vs_c2, method = "BH")
summary$padj_norm_vs_c3 <- p.adjust(summary$pval_norm_vs_c3, method = "BH")

write.csv(summary,"small_sum.csv")


######################    building simp heatmap data    ######################
######################       bottle neck on ram        ######################

### inner join sig methylation list  ###

sig_probes <- data.frame(beta_ex[,1])
colnames(sig_probes) <- c("probes")

notes2 <- notes1[order(notes1$cluster),]

write.csv(sig_probes, "sig_probes.csv")
write.csv(notes2, "notes2.csv")

###################   probe_to_gene.sh  ######################

manifest <- read_csv("headless.csv")

temp <- merge(x = sig_probes, y = manifest, by.x="probes", by.y="IlmnID", all = FALSE)
temp <- as.data.frame(temp$UCSC_RefGene_Name)
colnames(temp) <- c("genes")
genes <- temp %>%
  separate_rows(genes, sep=";")
conv_genes <- unique(genes)

write.csv(conv_genes,"conv_genes.csv")
