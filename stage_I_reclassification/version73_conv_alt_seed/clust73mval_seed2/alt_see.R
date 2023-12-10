library(tidyverse)

see <- read.csv("em_see.csv")
bval <- read_csv("beta_val_all.csv")
mval <- read_csv("mval_all.csv")

see2 <- subset(see, e0 == 0)
see2 <- see2[,-1]
see2$e0 <- 1

bval2 <- merge(x=see2, y=bval, by.x="master_list", by.y="file_names", all = FALSE)
bval2 <- bval2[,-c(2,3,4)]
rownames(bval2) <- bval2[,1]
colnames(bval2)[1] <- "file_names"

mval2 <- merge(x=see2, y=mval, by.x="master_list", by.y="...1", all = FALSE)
rownames(mval2) <- mval2[,1]
mval2 <- mval2[,-c(1,2)]

write.csv(see2, "em_see.csv")
write.csv(bval2, "beta_val_all.csv")
write.csv(mval2, "mval_all.csv")
