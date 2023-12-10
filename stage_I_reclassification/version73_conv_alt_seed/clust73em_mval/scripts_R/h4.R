
print ("loading significant methylation data")

temp <- read.csv("sig_meth_DE.csv")

cluster <- read.csv("clusters.csv")
temp2 <- t(cluster[,"temp"])
temp2 <- as.list(cbind("file_name", temp2))
colnames(temp) <- c(temp2)

temp2 <- temp[,1]
temp2 <- as.list(t(temp2))
rownames(temp) <- c(temp2)
temp <- temp[,-1]

notes2 <- read.csv("notes2.csv")
order <- notes2[,"File_Name"]

other <- temp[,order]

temp2 <- as.list(notes2$Sample_ID)
colnames(other) <- temp2
other2 <- log(other)

print ("building heatmap")

pdf("heatmap.pdf")
heatmap(as.matrix(other2), Colv = NA, scale="column", col = topo.colors(10))
dev.off()

########### building simplified heatmap ###########
print ("building simp heatmap...")

simp <- read.csv("simp_heat2.csv")
rownames(simp) <- as.list(t(simp[,1]))
simp <- simp[,-1]
simp <- log(simp)

pdf("simp_heatmap.pdf")
heatmap(as.matrix(simp), Colv = NA, scale="column", col = topo.colors(10))
dev.off()

pdf("simp_heatmap2.pdf")
heatmap(as.matrix(simp), Colv = NA, col = topo.colors(10))
dev.off()

legend <- c("ln0.1", "ln0.2", "ln0.3", "ln0.4", "ln0.5", "ln0.6", "ln0.7", "ln0.8", "ln0.9", "ln1.0")
pdf("legend.pdf")
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("center", legend=legend, fill= topo.colors(10))
dev.off()
