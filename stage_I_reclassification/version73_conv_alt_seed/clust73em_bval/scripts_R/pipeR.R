

# making sure that the bash scripts that are needed in this algorithm can be executed
command = "chmod +x ./pipe6.sh"
system(command)
command = "chmod +x ./bash_scripts/build_mx.sh"
system(command)
command = "chmod +x ./bash_scripts/new_emc.sh"
system(command)

# write up a function to do the within sum of squares needed for clustering
wssplot <- function(f1, nc=30, seed=1234){
  wss <- (nrow(f1)-1)*sum(apply(f1,2,var))
  for (i in 2:nc){
    set.seed(seed)
    wss[i] <- sum(kmeans(f1, centers=i)$withinss)}
  plot(1:nc, wss, type="b", xlab="Number of Clusters",
       ylab="Within groups sum of squares")
  wss
}

# new_best is just a binary counter to know if we need a replace old files
# it is also used to determine if the loop needs to keep going
new_best = 2

# num_sim is the number of simulation you would like to run in the loop
# counter counts down to know if the number of simulation all returned
# values worse than the initial EM
num_sim = 1000
counter = 1


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
    # em_see contains 2 rowns, a mater list of all probes, and the e0/best set of probes
    # first column contains all the probe ID
    # the second column contains only probe ID for the current best estimation
    # make a third column of 1 & 0, for the current best estimation
      # 1 --> the current best uses the probe in clustering
      # 2 --> the current best does NOT use the current probe in clustering
    em_see <- read.csv("em_see.csv", na.strings = c("", "NA"))
    em_see$e0 <- ifelse(is.na(em_see$Name),0,1)
    # we no longer need the second column, because the master list has our probe names
    em_see <- em_see[,-2]
    # we need to know the maximum number of probes possible, and the number of probes being used
    nprobes <- nrow(em_see)
    num_e0_probes <- sum(em_see$e0)
    
    # we set some value for the next estimation and clone our best estimate into a new column
    # we will remove up to half the current probes being used
    # we will add up to double probes being used
    # in most estimation this will leave us with more than our initial probes
    # but in some cases it will net us a loss of probes
    neg_probes = floor(runif(1, min=0, max=(num_e0_probes/2)))
    pos_probes = floor(runif(1, min=0, max=(min((num_e0_probes*5),nprobes))))
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
    write.csv(new,"probe_ex.csv")
    
    # remove the extra column that was create in the new probe estimate csv
    temp <- "sed -i 's/[^,]*,//' probe_ex.csv"
    system(temp)
    
    ### part 2 ###
    
    # pipe6 is a bash/sqlite3 script that will preform an inner join on the 
    # beta and m-value files. The beta and mval files contain a table each;
    # beta-val contains a table of i probes and j samples, each cell representing the beta score
    # where as mval contain cells representing the logit transformed beta value
    # the sqlite3 script will create new files with only the probes of interest
    print("step 2")
    command = "./pipe6.sh"
    system(command)
    
    # reading the m value file associated with only the probes of interest
    # pre-processing the files to remove any unwanted column taken up on 
    # the sql inner join and re-enter corrected files-names for collumn names
    # then transposing the matrix
    mx <- read.csv("mval_ex.csv")
    rownames(mx) = mx[,1]
    mx <- mx[,-c(1,2)]
    filenames <- read.csv("filenames.csv")
    colnames(mx) <- t(filenames)
    
    # set a seed and cluster the m-vales for the new estimation
    # we assume her 3 clusters, which might be subject of change/discussion later
    # create a csv that contains the file name and appropriate cluster
    set.seed(20)
    clusters <- as.data.frame(filenames)
    clusters$cluster <- kmeans(t(mx),3)$cluster
    write.csv(clusters,"clusters_ex.csv")
    
    # this bash/sqlite2 script creates a table with important information about
    # each sample, such as file name, case id, sample id, vital statue, 
    # and associate cluster
    # outputs a table for each cluster group using only the tumor sample of each
    # the bash script then uses shell command to create a list for each cluster,
    # the list has the number of samples, and the number of dead cases
    command = "./bash_scripts/build_mx.sh"
    system(command)
    
    ### part 3 ###
    
    # load the best current values of estimation
    print("step 3")
    na_best <- read.csv("na_best.csv")
    rownames(na_best) <- na_best[,1]
    na_best <- na_best[,-1]
    nb_best <- read.csv("nb_best.csv")
    rownames(nb_best) <- nb_best[,1]
    nb_best <- nb_best[,-1]
    
    # load up the list of sample;dead to create a table
    temp1 <- read.csv("cl1_nx.csv", header = FALSE)
    temp2 <- read.csv("cl2_nx.csv", header = FALSE)
    temp3 <- read.csv("cl3_nx.csv", header = FALSE)
    nx <- cbind(temp1, temp2, temp3)
    nx <- t(nx)
    rownames(nx) <- c("cl1", "cl2", "cl3")
    colnames(nx) <- c("samples", "dead", "live_score", "death_score")
    nx <- as.data.frame(nx)
    # find the ratio of dead patients / number of samples
    # this will be used to determine the better probes to use
    nx$ratio <- nx$dead/nx$samples
    nx$final_score <- nx$death_score - nx$live_score
    nx$rscore <- nx$ratio*(-1)*nx$final_score
    nxa <- nx[order(nx$final_score),]
    nxb <- nx[order(nx$rscore),]
    
    
    # complex if statement to decide what set of probes are best
    # if #1: worse cluster had a bestter death score, calc as n. deaths - (sum of: date of last contact of cencore samples / last recorded death day)
    # else if #2 worse cluster is the same, best cluster has a better ratio of death/sample
    # else if #3 same worse cluster, same best cluster, but less probes were used
    # else ... new estimation is worse than the current best estimation
    if ( nxa$final_score[3] > na_best$final_score[3] ) {
      # condition gives more accuracy results but not at the cost of finding less dead... nope.. that not how it worked out
      write.csv(nxa,"na_best.csv")
      write.csv(nxb,"nb_best.csv")
      print(paste(ex, " is best", sep=""))
      new_best = 1
    } else if ( round(nxa$final_score[3], digits = 5) == round(na_best$final_score[3], digits = 5) && nx$rscore[2] < nb_best$rscore[2] ) {
      # condition 1, the worse outcome is the same
      # but best outcome is better
      write.csv(nxa,"na_best.csv")
      write.csv(nxb,"nb_best.csv")
      print(paste(ex, " is best", sep=""))
      new_best = 1
    } else if ( round(nxa$final_score[3], digits = 5) == round(na_best$final_score[3], digits = 5) && round(nxb$rscore[2], digits = 5) == round(nb_best$rscore[2], digits = 5) && nrow(new) < num_e0_probes ) {
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
    print("step 4")
    if ( new_best == 1 ) {
      file.remove("current_probe.csv")
        file.copy(from = "probe_ex.csv",   # Copy files
                to = "current_probe.csv")

      file.remove("em_see.csv")
      command = "./bash_scripts/new_emc.sh"
      system(command)
    } else {
      counter <- counter - 1
      print(paste("end of trial ", y,"counter ", counter,  sep=""))
    }
    
    # clean up files to start a new simulation
    file.remove("probe_ex.csv")
    file.remove("cl1_nx.csv")
    file.remove("cl2_nx.csv")
    file.remove("cl3_nx.csv")
    file.remove("clusters_ex.csv")
    file.remove("beta_ex.csv")
    file.remove("mval_ex.csv")
    
    # mark down that we completed a siulated estimation
    y <- y+1
    
  }
}
