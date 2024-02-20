library(readr)
library(dplyr)

mval <- read.csv("mval_all.csv")
little <- mval[1:5,1:5]

mval2 <- mval
rownames(mval2) <- mval2[,1]
mval2 <- mval2[,-1]

filenames <- read.csv("filenames.csv")
colnames(mval2) <- filenames$names

filenames2 <- read.csv("filenames2.csv")

# Make sure that the 'names' column in df_secondary is character type
filenames2$CpG_filename <- as.character(filenames2$CpG_filename)

# Filter the columns of df_main
df_main_filtered <- mval2[, colnames(mval2) %in% filenames2$CpG_filename]

write.csv(df_main_filtered,"mval_all2.csv")
