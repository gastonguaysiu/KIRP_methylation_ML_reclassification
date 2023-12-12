library(tidyverse)

cl1_KIRP <- read.csv("cl1_GO.csv")
cl2_KIRP <- read.csv("cl2_GO.csv")
cl3_KIRP <- read.csv("cl3_GO.csv")
KIRP_GO <- read.csv("KIRP_GO.csv")

# Perform the inner join
cl1_KIRP <- inner_join(cl1_KIRP, KIRP_GO, by = c("X" = "X"))
cl2_KIRP <- inner_join(cl2_KIRP, KIRP_GO, by = c("X" = "X"))
cl3_KIRP <- inner_join(cl3_KIRP, KIRP_GO, by = c("X" = "X"))

cl1_KIRP2 <- cl1_KIRP[,1:11]
colnames(cl1_KIRP2) <- paste0("cl1_", colnames(KIRP_GO))

cl2_KIRP2 <- cl2_KIRP[,1:11]
colnames(cl2_KIRP2) <- paste0("cl2_", colnames(KIRP_GO))

cl3_KIRP2 <- cl3_KIRP[,1:11]
colnames(cl3_KIRP2) <- paste0("cl3_", colnames(KIRP_GO))

cl1cl2_cl3 <- inner_join(cl1_KIRP2, cl3_KIRP2, by = c("cl1_X" = "cl3_X"))

nx <- c(0, 0)
nx <- data.frame(nx)
rownames(nx) <- c("p_adj", "count")

nx[1,1] <- median(cl1_KIRP$p.adjust.x)
nx[2,1] <- median(cl1_KIRP$Count.x)
nx[1,2] <- median(cl1_KIRP$p.adjust.y)
nx[2,2] <- median(cl1_KIRP$Count.y)

nx[1,3] <- median(cl3_KIRP$p.adjust.x)
nx[2,3] <- median(cl3_KIRP$Count.x)
nx[1,4] <- median(cl3_KIRP$p.adjust.y)
nx[2,4] <- median(cl3_KIRP$Count.y)

nx[1,5] <- median(cl1cl2_cl3$cl1_p.adjust)
nx[2,5] <- median(cl1cl2_cl3$cl1_Count)
nx[1,6] <- median(cl1cl2_cl3$cl3_p.adjust)
nx[2,6] <- median(cl1cl2_cl3$cl3_Count)

colnames(nx) <- c("cl1", "KIRP", "cl3", "KIRP", "cl1_cl2_inner", "cl3_inner")

subset_df <- cl3_KIRP2[cl1_KIRP2$cl1_ONTOLOGY == "BP", ]
nrow(subset_df)

write.csv(nx,"GO_analysis.csv", row.names = FALSE)
