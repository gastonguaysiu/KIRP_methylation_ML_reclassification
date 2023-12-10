library("readr")


M0 <- read.csv("mval_M0.csv")
rownames(M0) = M0[,1]
M0 <- M0[,-c(1,2)]

clusters <- read.csv(("clusters.csv"))
colnames()


print ("cluster analysis...")

set.seed(20)

wssplot <- function(f1, nc=30, seed=1234){
  wss <- (nrow(f1)-1)*sum(apply(f1,2,var))
  for (i in 2:nc){
    set.seed(seed)
    wss[i] <- sum(kmeans(f1, centers=i)$withinss)}
  plot(1:nc, wss, type="b", xlab="Number of Clusters",
       ylab="Within groups sum of squares")
  wss
}

pdf("elbow1.pdf")
wssplot(M0, nc=10)
dev.off()

print ("writting clusters in promoter probes")
clusters2 <- as.data.frame(clusters[,1])
colnames(clusters2) <- "temp"
clusters2$cluster <- kmeans(t(M0),3)$cluster
write.csv(clusters2,"E0.csv", row.names = TRUE)
