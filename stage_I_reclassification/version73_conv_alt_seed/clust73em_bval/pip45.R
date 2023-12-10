
library(tidyverse)
library(plyr)

mval_all <- read_csv("mval_all.csv")
mval_ex <- merge(x = new, y = mval_all, by.x = "Name", by.y = "...1")

current_probe <- read.csv("current_probe.csv")
master <- read.csv("probe_e_all.csv")

em_see2 = em_see[,-2]
colnames(em_see2) = c("master_list", "e0")

write.csv(em_see2,"em_see2.csv")

em_see3 <- read.csv("em_see2.csv")
em_see3 <- em_see3[,-1]



ref <- read.csv("ref4.csv")
notes <- merge(x = ref, y = clusters, by.x = "File_Name", by.y = "names")
notes1 <- notes[,-c(4, 5, 9, 10, 11, 13, 14, 15)]

mxcl1 <- filter(notes1, cluster == 1, Sample_Type == "Primary Tumor")
mxcl1$live_score <- (mxcl1$lc_from_IPDD)/551.5
mxcl1$death_score <- (mxcl1$Days.to.death)/742.0

mxcl2 <- filter(notes1, cluster == 2, Sample_Type == "Primary Tumor")
mxcl2$live_score <- (mxcl2$lc_from_IPDD)/551.5
mxcl2$death_score <- (mxcl2$Days.to.death)/742.0

mxcl3 <- filter(notes1, cluster == 3, Sample_Type == "Primary Tumor")
mxcl3$live_score <- (mxcl3$lc_from_IPDD)/551.5
mxcl3$death_score <- (mxcl3$Days.to.death)/742.0

# initializing the data frame of new copy results
nxa <- c(0, 0, 0)
nxa <- data.frame(nxa)
rownames(nxa) <- c("cl1", "cl2", "cl3")
colnames(nxa) <- c("samples")
nxa$dead <- 0.01
nxa$live_weight <- 0
nxa$dead_weight = 0.01

nxa[1,1] <- nrow(mxcl1)
nxa[1,3] <- sum(mxcl1[, 'live_score'], na.rm = TRUE)
x <- length(which(mxcl1$"Vital.Status" == "Dead"))
if (x > 0 ) {
  nxa[1,2] <- x
  nxa[1,4] <- sum(mxcl1[, 'death_score'], na.rm = TRUE)
}

nxa[2,1] <- nrow(mxcl2)
nxa[2,3] <- sum(mxcl2[, 'live_score'], na.rm = TRUE)
x <- length(which(mxcl2$"Vital.Status" == "Dead"))
if (x > 0 ) {
  nxa[2,2] <- x
  nxa[2,4] <- sum(mxcl2[, 'death_score'], na.rm = TRUE)
}

nxa[3,1] <- nrow(mxcl3)
nxa[3,3] <- sum(mxcl3[, 'live_score'], na.rm = TRUE)
x <- length(which(mxcl3$"Vital.Status" == "Dead"))
if (x > 0 ) {
  nxa[3,2] <- x
  nxa[3,4] <- sum(mxcl3[, 'death_score'], na.rm = TRUE)
}



sum(mxcl1[, 'live_score'], na.rm = TRUE)
