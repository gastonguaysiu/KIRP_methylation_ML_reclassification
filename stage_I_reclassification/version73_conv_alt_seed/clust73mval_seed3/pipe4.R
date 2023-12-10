library(tidyverse)

# new_best is just a binary counter to know if we need a replace old files
# it is also used to determine if the loop needs to keep going
new_best = 2

# num_sim is the number of simulation you would like to run in the loop
# counter counts down to know if the number of simulation all returned
# values worse than the initial EM
num_sim = 10000
counter = 1

# let's speed things up by importing large tables only once
bval_all <- read_csv("beta_val_all.csv")
mval_all <- read_csv("mval_all.csv")

# if statement to get out of loop if we don't have improvements
while ( counter != 0 ) {
  y <- 1 
  counter <- num_sim
  print("outer while loop (line 37)")
  
  # the meat and bones of the EM algorith
  # y is a new counter just counting up until it hits the appropriate number of trials
  while (y < (num_sim + 1 )) {
    
    # ex is a string to represent the estimate trial
    ex <- paste("e",y,sep = "")
    print(paste("trial ", ex, " counter ", counter, sep=""))
    
    # em_see file contains a matrix of i rows each to represent a probe.
    # em_see contains 2 rows, a mater list of all probes, and the e0/best set of probes
    # first column contains all the probe ID
    # the second column contains 1s on probe ID for the current best estimation
    em_see <- read.csv("em_see.csv")
    em_see <- em_see[,-1]
    # we need to know the maximum number of probes possible, and the number of probes being used
    nprobes <- nrow(em_see)
    num_e0_probes <- sum(em_see$e0)
    
    # we set some value for the next estimation and clone our best estimate into a new column
    # we will remove up to half the current probes being used
    # we will add up to double probes being used
    # in most estimation this will leave us with more than our initial probes
    # but in some cases it will net us a loss of probes
    set.seed(floor(runif(1, min=0, max=(y*10))))
    neg_probes = floor(runif(1, min=0, max=(num_e0_probes/2)))
    pos_probes = floor(runif(1, min=0, max=(min((num_e0_probes/2),(nprobes-num_e0_probes)))))
    em_see$new <- em_see$e0
    
    # remove the probe in the clone estimate by setting 1 to 0 for probes being removed,
    # we choose a random number associated with the i-th probe and check if we are using the probe
    # this part might be able to be sub-optimal to run but at least it is definitely random
    while (pos_probes > 0) {
      x <- floor(runif(1, min=1, max=(nprobes+1)))
      if (em_see$new[x] == 0) { 
        em_see$new[x] <- 1
        pos_probes <- pos_probes - 1}
    }
    
    # add probes in the clone estimate by setting 0 to 1 for probes being added,
    # we choose a random number associated with the i-th probe and check if we not using the probe
    while (neg_probes > 0) {
      x <- floor(runif(1, min=1, max=(nprobes+1)))
      if (em_see$new[x] == 1) { 
        em_see$new[x] <- 0
        neg_probes <- neg_probes - 1}
    }
    
    # using the master list we create a new data frame that will contain the
    # names of the new estimated probes, which is then saved as a csv file
    # multi-step procedure to remove rows containing 0
    new <- cbind(em_see$master_list, em_see$new)
    new[new == 0] <- NA
    new <- na.omit(new)
    new <- as.data.frame(new[,1])
    colnames(new) <- "Name"
    
    ### part 2 ###
    
    # The beta and mval files contain a table each;
    # beta-val contains a table of i probes and j samples, each cell representing the beta score
    # where as mval contain cells representing the logit transformed beta value
    # the merge function acts like a inner join in this case
    beta_ex <- merge(x = new, y = bval_all, by.x = "Name", by.y = "...1")
    mx <- merge(x = new, y = mval_all, by.x = "Name", by.y = "...1")
    
    # pre-processing the files to remove any unwanted column taken up on 
    # the sql inner join and re-enter corrected files-names for column names
    # then transposing the matrix
    rownames(mx) = mx[,1]
    mx <- mx[,-1]
    filenames <- read.csv("filenames.csv")
    colnames(mx) <- t(filenames)
    
    # set a seed and cluster the m-vales for the new estimation
    # we assume her 3 clusters, which might be subject of change/discussion later
    # create a csv that contains the file name and appropriate cluster
    set.seed(3)
    clusters <- as.data.frame(filenames)
    clusters$cluster <- kmeans(t(mx),3,iter.max=10000)$cluster
    
    # create a table with important information about
    # each sample, such as file name, case id, sample id, vital statue, 
    # and associate cluster
    ref <- read.csv("ref4.csv")
    notes <- merge(x = ref, y = clusters, by.x = "File_Name", by.y = "names")
    notes1 <- notes[,-c(4, 5, 9, 10, 11, 13, 14, 15)]
    
    # create a set of tables for each cluster group using only the tumor sample of each
    mxcl1 <- filter(notes1, cluster == 1, Sample_Type == "Primary Tumor")
    mxcl1$live_score <- (mxcl1$lc_from_IPDD)/551.5
    mxcl1$death_score <- (mxcl1$Days.to.death)/742.0
    
    mxcl2 <- filter(notes1, cluster == 2, Sample_Type == "Primary Tumor")
    mxcl2$live_score <- (mxcl2$lc_from_IPDD)/551.5
    mxcl2$death_score <- (mxcl2$Days.to.death)/742.0
    
    mxcl3 <- filter(notes1, cluster == 3, Sample_Type == "Primary Tumor")
    mxcl3$live_score <- (mxcl3$lc_from_IPDD)/551.5
    mxcl3$death_score <- (mxcl3$Days.to.death)/742.0
    
    
    # create a minimal table to compare the results
    # initializing the data frame of new copy results
    nx <- c(0, 0, 0)
    nx <- data.frame(nx)
    rownames(nx) <- c("cl1", "cl2", "cl3")
    colnames(nx) <- c("samples")
    nx$dead <- 0.01
    nx$live_weight <- 0
    nx$dead_weight = 0.01
    
    nx[1,1] <- nrow(mxcl1)
    nx[1,3] <- sum(mxcl1[, 'live_score'], na.rm = TRUE)
    x <- length(which(mxcl1$"Vital.Status" == "Dead"))
    if (x > 0 ) {
      nx[1,2] <- x
      nx[1,4] <- sum(mxcl1[, 'death_score'], na.rm = TRUE)
    }
    
    nx[2,1] <- nrow(mxcl2)
    nx[2,3] <- sum(mxcl2[, 'live_score'], na.rm = TRUE)
    x <- length(which(mxcl2$"Vital.Status" == "Dead"))
    if (x > 0 ) {
      nx[2,2] <- x
      nx[2,4] <- sum(mxcl2[, 'death_score'], na.rm = TRUE)
    }
    
    nx[3,1] <- nrow(mxcl3)
    nx[3,3] <- sum(mxcl3[, 'live_score'], na.rm = TRUE)
    x <- length(which(mxcl3$"Vital.Status" == "Dead"))
    if (x > 0 ) {
      nx[3,2] <- x
      nx[3,4] <- sum(mxcl3[, 'death_score'], na.rm = TRUE)
    }
    
    ### part 3 ###

    # find the ratio of dead patients / number of samples
    # this will be used to determine the better probes to use
    nx$ratio <- nx$dead/nx$samples
    nx$final_score <- nx$dead - nx$live_weight
    nx$rscore <- (nx$ratio*nx$dead_weight)
    nxa <- nx[order(nx$final_score),]
    nxb <- nx[order(nx$rscore),]
    
    # load the best current values of estimation
    
    na_best <- read.csv("na_best.csv")
    rownames(na_best) <- na_best[,1]
    na_best <- na_best[,-1]
    nb_best <- read.csv("nb_best.csv")
    rownames(nb_best) <- nb_best[,1]
    nb_best <- nb_best[,-1]
    
    # complex if statement to decide what set of probes are best
    # if #1: worse cluster had a better death score, calc as n. deaths - (sum of: date of last contact of cencore samples / last recorded death day)
    # else if #2 worse cluster is the same, best cluster has a better overall survival (and accounts for at least 5% of samples => 311*0.05 = 15.55)
    # else if #3 same worse cluster, same best cluster, but less probes were used (but we want at least 1000 probes)
    # else ... new estimation is worse than the current best estimation
    if ( round(nxa$final_score[3], digits = 5) > round(na_best$final_score[3], digits = 5) ) {
      # condition 1, the worse outcome is the same
      # but best outcome is better
      write.csv(nxa,"na_best.csv")
      write.csv(nxb,"nb_best.csv")
      print(paste(ex, " is best", sep=""))
      new_best = 1
    } else if ( round(nxa$final_score[3], digits = 5) >= round(na_best$final_score[3], digits = 5) && round(nxb$rscore[1], digits = 7) < round(nb_best$rscore[1], digits = 7) ) {
      # best and worse is the same
      # but let's use less probes
      write.csv(nxa,"na_best.csv")
      write.csv(nxb,"nb_best.csv")
      print(paste(ex, " is best", sep=""))
      new_best = 1
    } else if ( round(nxa$final_score[3], digits = 5) >= round(na_best$final_score[3], digits = 5) && round(nxb$rscore[1], digits = 7) <= round(nb_best$rscore[1], digits = 7) && nrow(new) < num_e0_probes ) {
      # best and worse is the same
      # but let's use less probes
      write.csv(nxa,"na_best.csv")
      write.csv(nxb,"nb_best.csv")
      print(paste(ex, " is best", sep=""))
      new_best = 1
    } else {
      new_best = 0
    }
    
    ### part 4 ###
    
    # if the new estimation is better than the current best estimation,
      # replace the "best" files with the current estimate
    # else mark down the counter for not having found a better estimate
      # if all estimates in the simulation are are worse than the current best
      # loop ends, and we assume we have convergence
    if ( new_best == 1 ) {
      em_see2 = em_see[,-2]
      colnames(em_see2) = c("master_list", "e0")
      write.csv(em_see2,"em_see.csv")
    } else {
      counter <- counter - 1
    }
    
    # mark down that we completed a simulated estimation
    y <- y+1
    
  }
}
