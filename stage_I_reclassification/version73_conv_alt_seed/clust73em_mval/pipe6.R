library(tidyverse)

new_best <- 2
num_sim <- 10000
counter <- 1

bval_all <- read_csv("beta_val_all.csv")
mval_all <- read_csv("mval_all.csv")

while (counter != 0) {
  y <- 1
  counter <- num_sim

  while (y <= num_sim) {
    ex <- paste("e", y, sep = "")
    print(paste("trial ", ex, " counter ", counter, sep=""))

    em_see <- read.csv("em_see.csv") %>% select(-1)
    nprobes <- nrow(em_see)
    num_e0_probes <- sum(em_see$e0)

    set.seed(floor(runif(1, min=0, max=(y*10))))
    neg_probes <- floor(runif(1, min=0, max=(min(10 ,(num_e0_probes/2)))))
    pos_probes <- floor(runif(1, min=0, max=(min(10 ,(nprobes-num_e0_probes)))))
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
    
    set.seed(20)
    clusters <- as.data.frame(colnames(mx))
    colnames(clusters) <- "names"
    clusters$cluster <- kmeans(t(mx),3,iter.max=10000)$cluster

    ref <- read.csv("ref4.csv")
    notes <- merge(x = ref, y = clusters, by.x = "File_Name", by.y = "names")
    notes1 <- notes[,-c(4, 5, 9, 10, 11, 13, 14, 15)]

    mxcl_list <- lapply(1:3, function(x) {
      mxcl <- filter(notes1, cluster == x, Sample_Type == "Primary Tumor")
      mxcl$live_score <- (mxcl$lc_from_IPDD)/551.5
      mxcl$death_score <- (mxcl$Days.to.death)/742.0
      mxcl
    })

    nx <- data.frame(samples = c(0, 0, 0), dead = 0.01, live_weight = 0, dead_weight = 0.01)
    rownames(nx) <- c("cl1", "cl2", "cl3")

    for (i in 1:3) {
      mxcl <- mxcl_list[[i]]
      nx[i, "samples"] <- nrow(mxcl)
      nx[i, "live_weight"] <- sum(mxcl$live_score, na.rm = TRUE)
      x <- length(which(mxcl$Vital.Status == "Dead"))
      if (x > 0) {
        nx[i, "dead"] <- x
        nx[i, "dead_weight"] <- sum(mxcl$death_score, na.rm = TRUE)
      }
    }

    nx$ratio <- nx$dead/nx$samples
    nx$final_score <- nx$dead - nx$live_weight
    nx$rscore <- (nx$ratio*nx$dead_weight)
    nxa <- nx[order(nx$final_score),]
    nxb <- nx[order(nx$rscore),]

    na_best <- read.csv("na_best.csv")
    rownames(na_best) <- na_best[,1]
    na_best <- na_best[,-1]
    nb_best <- read.csv("nb_best.csv")
    rownames(nb_best) <- nb_best[,1]
    nb_best <- nb_best[,-1]

    if (round(nxa$final_score[3], digits = 5) > round(na_best$final_score[3], digits = 5)) {
      write.csv(nxa,"na_best.csv")
      write.csv(nxb,"nb_best.csv")
      print(paste0(ex, " is best, more accurate"))
      new_best = 1
    } else if (round(nxa$final_score[3], digits = 5) >= round(na_best$final_score[3], digits = 5) && round(nxb$rscore[1], digits = 9) < round(nb_best$rscore[1], digits = 9)) {
      write.csv(nxa,"na_best.csv")
      write.csv(nxb,"nb_best.csv")
      print(paste0(ex, " is best, more accurate"))
      new_best = 1
    } else if (round(nxa$final_score[3], digits = 5) >= round(na_best$final_score[3], digits = 5) && round(nxb$rscore[1], digits = 9) <= round(nb_best$rscore[1], digits = 9) && nrow(new) < num_e0_probes) {
      write.csv(nxa,"na_best.csv")
      write.csv(nxb,"nb_best.csv")
      print(paste0(ex, " is best, less probes"))
      new_best = 1
    } else {
      new_best = 0
    }

    if (new_best == 1) {
      em_see2 = em_see[,-2]
      colnames(em_see2) = c("master_list", "e0")
      write.csv(em_see2,"em_see.csv")
    } else {
      counter <- counter - 1
    }

    y <- y + 1
  }
}

print("second phase of EM algorithm, reintegrating probes")

new_best <- 2
num_sim <- 10000
counter <- 1

while (counter != 0) {
  y <- 1
  counter <- num_sim
  
  while (y <= num_sim) {
    ex <- paste("e", y, sep = "")
    print(paste("second phase, trial ", ex, " counter ", counter, sep=""))
    
    em_see <- read.csv("em_see.csv") %>% select(-1)
    nprobes <- nrow(em_see)
    num_e0_probes <- sum(em_see$e0)
    
    set.seed(floor(runif(1, min=0, max=(y*10))))
    neg_probes <- floor(runif(1, min=0, max=(min(10 ,(num_e0_probes/2)))))
    pos_probes <- floor(runif(1, min=0, max=min(10 ,(nprobes-num_e0_probes))))
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
    
    set.seed(20)
    clusters <- as.data.frame(colnames(mx))
    colnames(clusters) <- "names"
    clusters$cluster <- kmeans(t(mx),3,iter.max=10000)$cluster
    
    ref <- read.csv("ref4.csv")
    notes <- merge(x = ref, y = clusters, by.x = "File_Name", by.y = "names")
    notes1 <- notes[,-c(4, 5, 9, 10, 11, 13, 14, 15)]
    
    mxcl_list <- lapply(1:3, function(x) {
      mxcl <- filter(notes1, cluster == x, Sample_Type == "Primary Tumor")
      mxcl$live_score <- (mxcl$lc_from_IPDD)/551.5
      mxcl$death_score <- (mxcl$Days.to.death)/742.0
      mxcl
    })
    
    nx <- data.frame(samples = c(0, 0, 0), dead = 0.01, live_weight = 0, dead_weight = 0.01)
    rownames(nx) <- c("cl1", "cl2", "cl3")
    
    for (i in 1:3) {
      mxcl <- mxcl_list[[i]]
      nx[i, "samples"] <- nrow(mxcl)
      nx[i, "live_weight"] <- sum(mxcl$live_score, na.rm = TRUE)
      x <- length(which(mxcl$Vital.Status == "Dead"))
      if (x > 0) {
        nx[i, "dead"] <- x
        nx[i, "dead_weight"] <- sum(mxcl$death_score, na.rm = TRUE)
      }
    }
    
    nx$ratio <- nx$dead/nx$samples
    nx$final_score <- nx$dead - nx$live_weight
    nx$rscore <- (nx$ratio*nx$dead_weight)
    nxa <- nx[order(nx$final_score),]
    nxb <- nx[order(nx$rscore),]
    
    na_best <- read.csv("na_best.csv")
    rownames(na_best) <- na_best[,1]
    na_best <- na_best[,-1]
    nb_best <- read.csv("nb_best.csv")
    rownames(nb_best) <- nb_best[,1]
    nb_best <- nb_best[,-1]
    
    if (round(nxa$final_score[3], digits = 5) > round(na_best$final_score[3], digits = 5)) {
      write.csv(nxa,"na_best.csv")
      write.csv(nxb,"nb_best.csv")
      print(paste0(ex, " is best"))
      new_best = 1
    } else if (round(nxa$final_score[3], digits = 5) >= round(na_best$final_score[3], digits = 5) && round(nxb$rscore[1], digits = 9) < round(nb_best$rscore[1], digits = 9)) {
      write.csv(nxa,"na_best.csv")
      write.csv(nxb,"nb_best.csv")
      print(paste0(ex, " is best"))
      new_best = 1
    } else if (round(nxa$final_score[3], digits = 5) >= round(na_best$final_score[3], digits = 5) && round(nxb$rscore[1], digits = 9) <= round(nb_best$rscore[1], digits = 9) && nrow(new) > num_e0_probes) {
      write.csv(nxa,"na_best.csv")
      write.csv(nxb,"nb_best.csv")
      print(paste0(ex, " is best"))
      new_best = 1
    } else {
      new_best = 0
    }
    
    if (new_best == 1) {
      em_see2 = em_see[,-2]
      colnames(em_see2) = c("master_list", "e0")
      write.csv(em_see2,"em_see.csv")
    } else {
      counter <- counter - 1
    }
    
    y <- y + 1
  }
}
