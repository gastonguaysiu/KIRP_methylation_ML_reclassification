### just for testing on local
library("readr")
beta <- read_csv("beta_val_all.csv")

temp <- beta[,2]
temp <- as.list(t(temp))
beta <- beta[,-c(1,2)]
rownames(beta) <- c(temp)

cluster <- read.csv("clusters.csv")

temp <- as.list(cluster$temp)
colnames(beta) <- temp

beta_t <- t(beta)

temp <- as.data.frame(rownames(beta_t))
colnames(temp) <- "file_name"

beta_t <- cbind(temp, beta_t)

### building M0 data frame, this will be all the differentially methylated probe between Normal (living) subject and tumor (dead) subjects ###

print ("building a data frame for normal + each cluster")

normal <- read.csv("norm_L.csv")
normal <- merge(x=normal, y=beta_t, by.x="file_name", by.y="file_name")
temp <- as.list(t(normal[,1]))
rownames(normal) <- c(temp)
normal <- t(normal[,-1])

tumor <- read.csv("tumor_D.csv")
tumor  <- merge(x=tumor , y=beta_t, by.x="file_name", by.y="file_name")
temp <- as.list(t(tumor[,1]))
rownames(tumor) <- c(temp)
tumor <- t(tumor[,-1])


# create a function to be used in finding the standard deviation... fuck r base ##
rowVars <- function(x, na.rm=F) {
  # Vectorised version of variance filter
  rowSums((x - rowMeans(x, na.rm=na.rm))^2, na.rm=na.rm) / (ncol(x) - 1)
}

### finally getting those stats ###
print ("getting stats for norm + clusters...")

summary <- as.data.frame(rowMeans(normal))
colnames(summary) <- ("norm_avg")
summary$norm_std <- sqrt(rowVars(normal))
summary$tumor_avg <- rowMeans(tumor)
summary$tumor_std <- sqrt(rowVars(tumor))

### use Kolmogorov–Smirnov test of two samples (cluster vs norm) to get a p-value. KS tests shape of distribution unlike the t.test counterpart
# use Benjamini & Hochberg (1995) fdr, because stat quest youtube video explained this particular fdr

summary$pval <- sapply(1:nrow(normal), function(i) ks.test(as.vector(normal[i,]), as.vector(tumor[i,]))$p)
summary$padj <- p.adjust(summary$pval, method = "BH")

### create some tables to export to sqlite3, where we restrict to only significant probes
### i.e difference in beta expression of 25% and a p-adjusted of less than 0.01

M0 <- as.data.frame((pmax(summary$tumor_avg,summary$norm_avg))/(pmin(summary$tumor_avg,summary$norm_avg)))
colnames(M0) <- ("diff")
rownames(M0) <- rownames(summary)
M0$p_adj <- summary$padj
write.csv(M0,"M0.csv")
