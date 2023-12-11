
print ("loading significant methylation data")

########### --- ###########

temp <- read.csv("sig_meth_DE.csv", header = TRUE)

rownames(temp) <- temp[,1]
temp <- temp[,-1]
filenames <- read.csv("filenames.csv")
colnames(temp) <- t(filenames)
order_temp <- read.csv("notes2.csv")
order_temp <- order_temp[,-1]
order <- order_temp$File_Name

other <- temp[,order]

colnames(other) <- t(order_temp$Sample_ID)
other2 <- log(other)

########### dev. heatmap ###########

print ("building simp heatmap...")

simp <- read.csv("simp_heat2.csv")
rownames(simp) <- as.list(t(simp[,1]))
simp <- simp[,-1]
simp <- log(simp)

simp2 <- simp
colnames(simp2) <- c("t_cl1", "t_cl2", "n_cl3", "t_cl3")

########### dev. ComplexHeatmap ###########

library(ComplexHeatmap)

mat = as.matrix(simp2[,c(1,2,4,3)])
ha = HeatmapAnnotation(grouping = c(1, 2, 3, 3))
Heatmap(mat, name = "m_value", show_row_names = FALSE, col = topo.colors(10), top_annotation = ha)
Heatmap(mat, name = "m_value", show_row_names = FALSE, cluster_columns = FALSE, col = topo.colors(10), top_annotation = ha)

simp2 <- simp2[,c(1,2,4,3)]
heatmap(as.matrix(simp2), Colv = NA, col = topo.colors(10),
        main = "Average ln(B) value methylation score",
        sub ="n = normal tissue sample, t = Primary tumor sample, Beta values averaged together then transformed")


mat = as.matrix(other2)
an <- order_temp[,9]
ha = HeatmapAnnotation(grouping = an)
Heatmap(mat, name = "m_value", show_row_names = FALSE, show_column_names = FALSE, col = topo.colors(10), top_annotation = ha)
Heatmap(mat, name = "m_value", show_row_names = FALSE, show_column_names = FALSE, cluster_columns = FALSE, col = topo.colors(10), top_annotation = ha)



