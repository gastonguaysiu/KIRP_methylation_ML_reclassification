em_see <- read.csv("em_see.csv", na.strings = c("", "NA"))

em_see$e0 <- ifelse(is.na(em_see$Name),0,1)
em_see <- em_see[,-2]

nprobes <- nrow(em_see)
num_e0_probes <- sum(em_see$e0)

neg_probes = floor(runif(1, min=0, max=(num_e0_probes/2)))
pos_probes = floor(runif(1, min=0, max=(num_e0_probes)))

y <- 1 
while (y < 2) {
  ex <- paste("e",y,sep = "")
  em_see$new <- em_see$e0
  
  neg_probes = floor(runif(1, min=0, max=(num_e0_probes/2)))
  pos_probes = floor(runif(1, min=0, max=(num_e0_probes)))
  em_see$new <- em_see$e0
  
  while (pos_probes > 0) {
    x <- floor(runif(1, min=1, max=(nprobes+1)))
    if (em_see$new[x] == 0) { 
      em_see$new[x] <- 1
      pos_probes <- pos_probes - 1}
  }
  
  while (neg_probes > 0) {
    x <- floor(runif(1, min=1, max=(nprobes+1)))
    if (em_see$new[x] == 1) { 
      em_see$new[x] <- 0
      neg_probes <- neg_probes - 1}
  }
  
  new <- cbind(em_see$master_list, em_see$new)
  new[new == 0] <- NA
  new <- na.omit(new)
  new <- as.data.frame(new[,1])
  colnames(new) <- "Name"
  write.csv(new,file=paste("probe_", ex, ".csv", sep=""))
  
  temp <- paste("sed -i 's/[^,]*,//' probe_" , ex, ".csv", sep="")
  system(temp)
  
  names(em_see)[names(em_see) == 'new'] <- ex
  y <- y+1
}

# temp <- paste("sed -i 's/[^,]*,//' probe_" , ex, ".csv", sep="")
# cat(temp)
