
y <- 1 
while (y < 2) {
  
  ex <- paste("e",y,sep = "")
  
  n_best <- read.csv("n_best.csv")
  rownames(n_best) <- n_best[,1]
  n_best <- n_best[,-1]
  
  temp1 <- read.csv(file=paste("cl1_n", y, ".csv", sep=""), header = FALSE)
  temp2 <- read.csv(file=paste("cl2_n", y, ".csv", sep=""), header = FALSE)
  temp3 <- read.csv(file=paste("cl3_n", y, ".csv", sep=""), header = FALSE)
  nx <- cbind(temp1, temp2, temp3)
  nx <- t(nx)
  rownames(nx) <- c("cl1", "cl2", "cl3")
  colnames(nx) <- c("samples", "dead")
  nx <- as.data.frame(nx)
  nx$ratio <- nx$dead/nx$samples
  
  nx <- nx[order(nx$ratio),]

  if ( nx$ratio[3] > n_best$ratio[3] && nx$dead[3] >= n_best$dead[3] ) {
    # condition gives better equal or more accuracy but not at the cost of finding less dead
    n_best <- nx
    write.csv(n_best,"n_best.csv")
    print(paste(ex, " is best", sep=""))
    new_best = 1
  } else if ( nx$ratio[3] == n_best$ratio[3] && nx$dead[3] >= n_best$dead[3] && nx$ratio[2] >= n_best$ratio[2]) {
    # condition 1 & 2 means the worst outcome group is about the saem, 
    # but the middle outcome is equally or more accurate
    n_best <- nx
    write.csv(n_best,"n_best.csv")
    print(paste(ex, " is best", sep=""))
    new_best = 1
  } else {
    new_best = 0
  }
  
  if ( new_best = 1 ) {
    file.remove("current_probe.csv")
    file.copy(from = paste(file = "probe_e", y, ".csv", sep=""),   # Copy files
              to = "current_probe.csv")
  }
file.remove(paste("probe_e", y, ".csv", sep=""))
file.remove(paste("cl1_n", y, ".csv", sep=""))
file.remove(paste("cl2_n", y, ".csv", sep=""))
file.remove(paste("cl3_n", y, ".csv", sep=""))
file.remove(paste("clusters_e", y, ".csv", sep=""))
file.remove(paste("beta_e", y, ".csv", sep=""))
file.remove(paste("mval_e", y, ".csv", sep=""))
  }

  y <- y+1
}
