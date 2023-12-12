librart(tidyverse)
library(stats)

fisher <- read.csv("fisher_stats.csv")
rownames(fisher) <- fisher[,1]
fisher <- fisher[,-1]
fisher <- t(fisher)


fisher_cl1 <- fisher.test(fisher[,c(1,2)], simulate.p.value = TRUE, B = 10000)$p.value
fisher_cl2 <- fisher.test(fisher[,c(3,4)], simulate.p.value = TRUE, B = 10000)$p.value
fisher_cl3 <- fisher.test(fisher[,c(5,6)], simulate.p.value = TRUE, B = 10000)$p.value

df <- data.frame(0, 0, 0, 0, 0, 0)
df[1,] <- fisher_cl1 
df[,c(3,4)] <- fisher_cl2 
df[,c(5,6)] <- fisher_cl3
rownames(df) <- c("fisher_pval")
colnames(df) <- colnames(fisher)

fisher <- t(rbind(fisher, df))
write.csv(fisher,"fisher_val.csv")
