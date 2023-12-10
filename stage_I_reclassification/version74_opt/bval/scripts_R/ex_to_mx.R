

command = "chmod +x ./pipe6.sh"
system(command)

command = "./pipe6.sh"
system(command)

wssplot <- function(f1, nc=30, seed=1234){
  wss <- (nrow(f1)-1)*sum(apply(f1,2,var))
  for (i in 2:nc){
    set.seed(seed)
    wss[i] <- sum(kmeans(f1, centers=i)$withinss)}
  plot(1:nc, wss, type="b", xlab="Number of Clusters",
       ylab="Within groups sum of squares")
  wss
}

y <- 1 
while (y < 2) {
  mval_mx <- paste("mval_e",y,".csv", sep = "")
  
  mx <- read.csv(mval_mx)
  rownames(mx) = mx[,1]
  mx <- mx[,-c(1,2)]
  filenames <- read.csv("filenames.csv")
  colnames(mx) <- t(filenames)
  
  set.seed(20)
  
  clusters <- as.data.frame(filenames)
  clusters$cluster <- kmeans(t(mx),3)$cluster
  write.csv(clusters,file=paste("clusters_e", y, ".csv", sep=""))

  y <- y+1
}


command = "chmod +x ./bash_scripts/build_mx.sh"
system(command)

command = "./bash_scripts/build_mx.sh"
system(command)
