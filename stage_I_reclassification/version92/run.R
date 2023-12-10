library(tidyverse)

# write up a function to do the within sum of squares needed for clustering
wssplot <- function(f1, nc=30, seed=1234){
  wss <- (nrow(f1)-1)*sum(apply(f1,2,var))
  for (i in 2:nc){
    set.seed(seed)
    wss[i] <- sum(kmeans(f1, centers=i)$withinss)}
  plot(1:nc, wss, type="b", xlab="Number of Clusters",
       ylab="Within groups sum of squares")
  wss
}

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

# set a seed and cluster the m-vales for the new estimation
# we assume her 3 clusters, which might be subject of change/discussion later
# create a csv that contains the file name and appropriate cluster
pdf("elbow1.pdf")
wssplot(mx, nc=10)
dev.off()

set.seed(20)
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

# write.csv(mx, "mval_ex.csv")
# write.csv(bx, "bval_ex.csv")
# write.csv(clusters,"clusters.csv")
#
# write.csv(cluster1, "cluster1_info.csv")
# write.csv(cluster2, "cluster2_info.csv")
# write.csv(cluster3, "cluster3_info.csv")
# write.csv(normal, "normal_info.csv")
# write.csv(normal1, "normal1_info.csv")
# write.csv(normal2, "normal2_info.csv")
# write.csv(normal3, "normal3_info.csv")

###################   stats.R   ######################

print ("pre-processing and loading up beta_current...")

beta_current <- read.csv("beta_promo.csv")

temp <- beta_current[,3]
temp <- as.list(t(temp))
beta_current <- beta_current[,-c(1,2,3)]
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

write.csv(summary,"stat_summary.csv")

### create some tables to export to sqlite3, where we restrict to only significant probes
### i.e difference in beta expression of 33% and a p-adjusted of less than 0.01

probes_cl1 <- as.data.frame((summary$c1_avg)/(summary$norm_avg))
colnames(probes_cl1) <- ("diff")
rownames(probes_cl1) <- rownames(summary)
probes_cl1$p_adj <- summary$padj_norm_vs_c1

probes_cl2 <- as.data.frame((summary$c2_avg)/(summary$norm_avg))
colnames(probes_cl2) <- ("diff")
rownames(probes_cl2) <- rownames(summary)
probes_cl2$p_adj <- summary$padj_norm_vs_c2

probes_cl3 <- as.data.frame((summary$c3_avg)/(summary$norm_avg))
colnames(probes_cl3) <- ("diff")
rownames(probes_cl3) <- rownames(summary)
probes_cl3$p_adj <- summary$padj_norm_vs_c3

# write.csv(probes_cl1,"probe_cl1.csv")
# write.csv(probes_cl2,"probe_cl2.csv")
# write.csv(probes_cl3,"probe_cl3.csv")

### separating info into data frames to later in a simplified heatmap ###

print ("building a data frame for normal an/or cluster for simp")
print ("please edit this section before running, not all clusters will have normal samples")

# norm_cl1 <- read.csv("norm_cl1.csv")
# norm_cl1 <- merge(x=norm_cl1, y=beta_t, by.x="File_Name", by.y="file_name")
# temp <- as.list(t(norm_cl1[,1]))
# rownames(norm_cl1) <- c(temp)
# norm_cl1 <- t(norm_cl1[,-1])

# norm_cl2 <- read.csv("norm_cl2.csv")
# norm_cl2 <- merge(x=norm_cl2, y=beta_t, by.x="File_Name", by.y="file_name")
# temp <- as.list(t(norm_cl2[,1]))
# rownames(norm_cl2) <- c(temp)
# norm_cl2 <- t(norm_cl2[,-1])

norm_cl3 <- as.data.frame(normal3$File_Name)
colnames(norm_cl3) <- c("File_Name")
norm_cl3 <- merge(x=norm_cl3, y=beta_t, by.x="File_Name", by.y="file_name", all = FALSE)
temp <- as.list(t(norm_cl3[,1]))
rownames(norm_cl3) <- c(temp)
norm_cl3 <- t(norm_cl3[,-1])

### finally getting averages to build a simplified heatmap ###
print ("getting averages for norm and/or clusters...")

# simp_heat <- as.data.frame(rowMeans(norm_cl1))
# colnames(simp_heat) <- ("norm_cl1_avg")
simp_heat <- as.data.frame(rowMeans(cluster1B))
colnames(simp_heat) <- ("tumor_cl1_avg")
# simp_heat$norm_cl2_avg <- rowMeans(norm_cl2)
simp_heat$tumor_cl2_avg <- rowMeans(cluster2B)
simp_heat$norm_cl3_avg <- rowMeans(norm_cl3)
simp_heat$tumor_cl3_avg <- rowMeans(cluster3B)

###################   sig_probe.sh   ######################

cl1B <- probes_cl1[(probes_cl1$diff > 1.33 | probes_cl1$diff < 0.67) & probes_cl1$p_adj < 0.01,]
cl2B <- probes_cl2[(probes_cl2$diff > 1.33 | probes_cl2$diff < 0.67) & probes_cl2$p_adj < 0.01,]
cl3B <- probes_cl3[(probes_cl3$diff > 1.33 | probes_cl3$diff < 0.67) & probes_cl3$p_adj < 0.01,]

cl1b <- as.data.frame(rownames(cl1B))
cl2b <- as.data.frame(rownames(cl2B))
cl3b <- as.data.frame(rownames(cl3B))

colnames(cl1b) <- c("probes")
colnames(cl2b) <- c("probes")
colnames(cl3b) <- c("probes")

sig_probes <- rbind(cl1b, cl2b, cl3b)
sig_probes <- as.data.frame(sig_probes[!duplicated(sig_probes), ])
colnames(sig_probes) <- c("probes")

######################    building simp heatmap data    ######################
######################       bottle neck on ram        ######################

### inner join sig methylation list  ###

beta_promo <- read_csv("beta_promo.csv")
beta_promo <- beta_promo[,-c(2,3)]
sig_meth_DE <- merge(x=sig_probes, y=beta_promo, by.x="probes", by.y="Name", all = FALSE)
rownames(sig_meth_DE) <- sig_meth_DE[,1]
sig_meth_DE <- sig_meth_DE[,-1]
colnames(sig_meth_DE) <- t(filenames)

simp_heat$temp <- rownames(simp_heat)
simp_heat2 <- merge(x=sig_probes, y=simp_heat, by.x="probes", by.y="temp", all = FALSE)
rownames(simp_heat2) <- simp_heat2[,1]
simp_heat2 <- simp_heat2[,-1]

notes2 <- notes1[order(notes1$cluster),]

write.csv(sig_meth_DE, "sig_meth_DE.csv")
write.csv(notes2, "notes2.csv")
write.csv(simp_heat2,"simp_heat2.csv")

########### dev. simp heatmaps ###########
#
# print ("building simp heatmap...")
#
# simp <- log(simp_heat2)
# colnames(simp) <- c("t_cl1", "t_cl2", "n_cl3", "t_cl3")
#
# library(ComplexHeatmap)
#
# mat = as.matrix(simp[,c(1,2,4,3)])
# ha = HeatmapAnnotation(grouping = c(1, 2, 3, 3))
# Heatmap(mat, name = "m_value", show_row_names = FALSE, col = topo.colors(10), top_annotation = ha)
# Heatmap(mat, name = "m_value", show_row_names = FALSE, cluster_columns = FALSE, col = topo.colors(10), top_annotation = ha)
#
# simp2 <- simp[,c(1,2,4,3)]
# heatmap(as.matrix(simp2), Colv = NA, col = topo.colors(10),
#         main = "Average ln(B) value methylation score",
#         sub ="n = normal tissue sample, t = Primary tumor sample, Beta values averaged together then transformed")
#

###################   probe_to_gene.sh  ######################

manifest <- read_csv("headless.csv")

cl_temp <- merge(x = cl1b, y = manifest, by.x="probes", by.y="IlmnID", all = FALSE)
cl_temp <- as.data.frame(cl_temp$UCSC_RefGene_Name)
colnames(cl_temp) <- c("genes")
cl_genes <- cl_temp %>%               
  separate_rows(genes, sep=";")
cl1_dm_genes <- unique(cl_genes)

cl_temp <- merge(x = cl2b, y = manifest, by.x="probes", by.y="IlmnID", all = FALSE)
cl_temp <- as.data.frame(cl_temp$UCSC_RefGene_Name)
colnames(cl_temp) <- c("genes")
cl_genes <- cl_temp %>%               
  separate_rows(genes, sep=";")
cl2_dm_genes <- unique(cl_genes)

cl_temp <- merge(x = cl3b, y = manifest, by.x="probes", by.y="IlmnID", all = FALSE)
cl_temp <- as.data.frame(cl_temp$UCSC_RefGene_Name)
colnames(cl_temp) <- c("genes")
cl_genes <- cl_temp %>%               
  separate_rows(genes, sep=";")
cl3_dm_genes <- unique(cl_genes)

###################   hyper hypo analysis  ######################

hyper_cl1B <- probes_cl1[probes_cl1$diff > 1.33 & probes_cl1$p_adj < 0.01,]
hypo_cl1B <- probes_cl1[probes_cl1$diff < 0.67 & probes_cl1$p_adj < 0.01,]
hyper_cl2B <- probes_cl2[probes_cl2$diff > 1.33 & probes_cl2$p_adj < 0.01,]
hypo_cl2B <- probes_cl2[probes_cl2$diff < 0.67 & probes_cl2$p_adj < 0.01,]
hyper_cl3B <- probes_cl3[probes_cl3$diff > 1.33 & probes_cl3$p_adj < 0.01,]
hypo_cl3B <- probes_cl3[probes_cl3$diff < 0.67 & probes_cl3$p_adj < 0.01,]

hyper_cl1b <- as.data.frame(rownames(hyper_cl1B))
hypo_cl1b <- as.data.frame(rownames(hypo_cl1B))
hyper_cl2b <- as.data.frame(rownames(hyper_cl2B))
hypo_cl2b <- as.data.frame(rownames(hypo_cl2B))
hyper_cl3b <- as.data.frame(rownames(hyper_cl3B))
hypo_cl3b <- as.data.frame(rownames(hypo_cl3B))

colnames(hyper_cl1b) <- c("probes")
colnames(hypo_cl1b) <- c("probes")
colnames(hyper_cl2b) <- c("probes")
colnames(hypo_cl2b) <- c("probes")
colnames(hyper_cl3b) <- c("probes")
colnames(hypo_cl3b) <- c("probes")

###################   probe_to_gene.sh  ######################

cl_temp <- merge(x = hyper_cl1b, y = manifest, by.x="probes", by.y="IlmnID", all = FALSE)
cl_temp <- as.data.frame(cl_temp$UCSC_RefGene_Name)
colnames(cl_temp) <- c("genes")
cl_genes <- cl_temp %>%
  separate_rows(genes, sep=";")
hyper_cl1_dm_genes <- unique(cl_genes)

cl_temp <- merge(x = hypo_cl1b, y = manifest, by.x="probes", by.y="IlmnID", all = FALSE)
cl_temp <- as.data.frame(cl_temp$UCSC_RefGene_Name)
colnames(cl_temp) <- c("genes")
cl_genes <- cl_temp %>%
  separate_rows(genes, sep=";")
hypo_cl1_dm_genes <- unique(cl_genes)

hyper_hypo_overlap_cl1 <- merge(x = hyper_cl1_dm_genes, y = hypo_cl1_dm_genes, by.x="genes", by.y="genes", all = FALSE)

cl_temp <- merge(x = hyper_cl2b, y = manifest, by.x="probes", by.y="IlmnID", all = FALSE)
cl_temp <- as.data.frame(cl_temp$UCSC_RefGene_Name)
colnames(cl_temp) <- c("genes")
cl_genes <- cl_temp %>%
  separate_rows(genes, sep=";")
hyper_cl2_dm_genes <- unique(cl_genes)

cl_temp <- merge(x = hypo_cl2b, y = manifest, by.x="probes", by.y="IlmnID", all = FALSE)
cl_temp <- as.data.frame(cl_temp$UCSC_RefGene_Name)
colnames(cl_temp) <- c("genes")
cl_genes <- cl_temp %>%
  separate_rows(genes, sep=";")
hypo_cl2_dm_genes <- unique(cl_genes)

hyper_hypo_overlap_cl2 <- merge(x = hyper_cl2_dm_genes, y = hypo_cl2_dm_genes, by.x="genes", by.y="genes", all = FALSE)

cl_temp <- merge(x = hyper_cl3b, y = manifest, by.x="probes", by.y="IlmnID", all = FALSE)
cl_temp <- as.data.frame(cl_temp$UCSC_RefGene_Name)
colnames(cl_temp) <- c("genes")
cl_genes <- cl_temp %>%
  separate_rows(genes, sep=";")
hyper_cl3_dm_genes <- unique(cl_genes)

cl_temp <- merge(x = hypo_cl3b, y = manifest, by.x="probes", by.y="IlmnID", all = FALSE)
cl_temp <- as.data.frame(cl_temp$UCSC_RefGene_Name)
colnames(cl_temp) <- c("genes")
cl_genes <- cl_temp %>%
  separate_rows(genes, sep=";")
hypo_cl3_dm_genes <- unique(cl_genes)

hyper_hypo_overlap_cl3 <- merge(x = hyper_cl3_dm_genes, y = hypo_cl3_dm_genes, by.x="genes", by.y="genes", all = FALSE)

######  make some files for GO enrichment analysis  ######

write.csv(cl1_dm_genes,"cl1_dm_genes.csv")
write.csv(hyper_cl1_dm_genes,"cl1_dm_genes_hyper.csv")
write.csv(hypo_cl1_dm_genes,"cl1_dm_genes_hypo.csv")
write.csv(cl2_dm_genes,"cl2_dm_genes.csv")
write.csv(hyper_cl2_dm_genes,"cl2_dm_genes_hyper.csv")
write.csv(hypo_cl2_dm_genes,"cl2_dm_genes_hypo.csv")
write.csv(cl3_dm_genes,"cl3_dm_genes.csv")
write.csv(hyper_cl3_dm_genes,"cl3_dm_genes_hyper.csv")
write.csv(hypo_cl3_dm_genes,"cl3_dm_genes_hypo.csv")

# lazy bash column correction
command = "sed -i 's/[^,]*,//' cl*"
system(command)
### sed 's/"//g' cl1_dm_genes.csv > cl1_dm_genes2.csv
