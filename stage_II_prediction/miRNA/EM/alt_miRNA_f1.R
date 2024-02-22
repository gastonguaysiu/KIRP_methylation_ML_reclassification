library(tidyverse)

filenames <- read.csv("filename_trainning.csv")

# Create folders conv1 to conv10 if they don't exist
for (i in 1:10) {
  
  miRNA <- read_csv("trainingData.csv")
  
  em_see <- read.csv("em_see.csv")[,-1]
  nprobes <- nrow(em_see)
  num_e0_probes <- sum(em_see$e0)
  
  nprobes <- nrow(em_see)
  num_e0_probes <- sum(em_see$e0)
  em_see$new <- em_see$e0
  
  new <- em_see %>%
    select(master_list, new) %>%
    mutate(new = if_else(new == 0, NA_real_, new)) %>%
    na.omit() %>%
    select(-new) %>%
    setNames("Name")
  
  mx <- merge(x = new, y = miRNA, by.x = "Name", by.y = "...1") %>%
    select(-1)
  
  set.seed(42)
  clusters <- tibble(sample = colnames(mx), 
                     cluster = kmeans(t(mx), 3, iter.max = 10000)$cluster)
  
  ref <- read.csv("gdc_join_clinical.csv")
  notes <- merge(x = ref, y = clusters, by.x = "miRNA_filename", by.y = "sample") %>%
    select(-c(5,6,15:28))
  
  calculate_scores <- function(df, cluster_id) {
    df <- df %>%
      filter(cluster == cluster_id, Sample_Type == "Primary Tumor") %>%
      mutate(live_score = as.numeric(days_to_last_follow_up) / 624,
             death_score = as.numeric(days_to_death) / 761)
  }
  
  
  mxcls <- lapply(1:3, calculate_scores, df = notes)
  
  nx <- data.frame(samples = vapply(mxcls, nrow, numeric(1)),
                   dead = vapply(mxcls, function(df) sum(df$vital_status == "Dead"), numeric(1)),
                   live_weight = vapply(mxcls, function(df) sum(df$live_score, na.rm = TRUE), numeric(1)),
                   dead_weight = 0.01)
  
  total_dead <- sum(nx$dead)
  
  nx <- nx %>%
    rownames_to_column("cluster") %>%
    mutate(dead = if_else(dead == 0, 0.01, dead),
           dead_weight = if_else(dead > 0, if_else(vapply(mxcls, function(df) sum(df$death_score, na.rm = TRUE), numeric(1)) == 0, 0.01, vapply(mxcls, function(df) sum(df$death_score, na.rm = TRUE), numeric(1))), dead_weight),
           ratio = dead / samples,
           precision = dead / (dead + live_weight),
           recall = dead / total_dead,
           final_score = if_else((precision + recall) == 0, 0, 2 * (precision * recall) / (precision + recall)),
           rscore = ratio * dead_weight)
  
  nxa <- nx %>% arrange(final_score)
  nxb <- nx %>% arrange(rscore)
  write.csv(nxa, "na_best.csv")
  write.csv(nxb, "nb_best.csv")
  
#########################################################################
  
  new_best <- 2
  num_sim <- 10000
  counter <- num_sim
  y <- 1
  
  miRNA <- read_csv("trainingData.csv")
  
  while ( counter != 0 ) {
    
    if ( y > 10000 ) {
      y <- 1
      counter <- num_sim
    }
    
    print(paste("trial e", y, " counter ", counter, sep = ""))
    
    em_see <- read.csv("em_see.csv")[,-1]
    nprobes <- nrow(em_see)
    num_e0_probes <- sum(em_see$e0)
    
    nprobes <- nrow(em_see)
    num_e0_probes <- sum(em_see$e0)
    
    set.seed(floor(runif(1, min=0, max=(y*10))))
    neg_probes = floor(runif(1, min=0, max=(num_e0_probes/2)))
    pos_probes = floor(runif(1, min=0, max=(min((num_e0_probes/2),(nprobes-num_e0_probes)))))
    em_see$new <- em_see$e0
    
    # find indices of 0s and 1s
    zero_indices <- which(em_see$new == 0)
    one_indices <- which(em_see$new == 1)
    
    # sample and update for adding probes
    if (length(zero_indices) >= pos_probes && pos_probes > 0) {
      indices_to_update <- sample(zero_indices, pos_probes)
      em_see$new[indices_to_update] <- 1
    }
    
    # sample and update for removing probes
    if (length(one_indices) >= neg_probes && neg_probes > 0) {
      indices_to_update <- sample(one_indices, neg_probes)
      em_see$new[indices_to_update] <- 0
    }
    
    new <- em_see %>%
      select(master_list, new) %>%
      mutate(new = if_else(new == 0, NA_real_, new)) %>%
      na.omit() %>%
      select(-new) %>%
      setNames("Name")
    
    mx <- merge(x = new, y = miRNA, by.x = "Name", by.y = "...1") %>%
      select(-1)
    
    set.seed(42)
    clusters <- tibble(sample = colnames(mx), 
                       cluster = kmeans(t(mx), 3, iter.max = 10000)$cluster)
    
    ref <- read.csv("gdc_join_clinical.csv")
    notes <- merge(x = ref, y = clusters, by.x = "miRNA_filename", by.y = "sample") %>%
      select(-c(5,6,15:28))
    
    calculate_scores <- function(df, cluster_id) {
      df <- df %>%
        filter(cluster == cluster_id, Sample_Type == "Primary Tumor") %>%
        mutate(live_score = as.numeric(days_to_last_follow_up) / 624,
               death_score = as.numeric(days_to_death) / 761)
    }
    
    
    mxcls <- lapply(1:3, calculate_scores, df = notes)
    
    nx <- data.frame(samples = vapply(mxcls, nrow, numeric(1)),
                     dead = vapply(mxcls, function(df) sum(df$vital_status == "Dead"), numeric(1)),
                     live_weight = vapply(mxcls, function(df) sum(df$live_score, na.rm = TRUE), numeric(1)),
                     dead_weight = 0.01)
    
    total_dead <- sum(nx$dead)
    
    nx <- nx %>%
      rownames_to_column("cluster") %>%
      mutate(dead = if_else(dead == 0, 0.01, dead),
             dead_weight = if_else(dead > 0, if_else(vapply(mxcls, function(df) sum(df$death_score, na.rm = TRUE), numeric(1)) == 0, 0.01, vapply(mxcls, function(df) sum(df$death_score, na.rm = TRUE), numeric(1))), dead_weight),
             ratio = dead / samples,
             precision = dead / (dead + live_weight),
             recall = dead / total_dead,
             final_score = if_else((precision + recall) == 0, 0, 2 * (precision * recall) / (precision + recall)),
             rscore = ratio * dead_weight)
    
    nxa <- nx %>% arrange(final_score)
    nxb <- nx %>% arrange(rscore)
    
    na_best <- read.csv("na_best.csv") %>% select(-1)
    nb_best <- read.csv("nb_best.csv") %>% select(-1)
    
    if (round(nxa$final_score[3], 5) > round(na_best$final_score[3], 5) |
        (round(nxa$final_score[3], 5) >= round(na_best$final_score[3], 5) & round(nxb$rscore[1], 7) < round(nb_best$rscore[1], 7)) |
        (round(nxa$final_score[3], 5) >= round(na_best$final_score[3], 5) & round(nxb$rscore[1], 7) <= round(nb_best$rscore[1], 7) & nrow(new) < num_e0_probes)) {
      write.csv(nxa, "na_best.csv")
      write.csv(nxb, "nb_best.csv")
      
      print(paste("e", y, " is best", sep = ""))
      file.remove("em_see.csv")
      em_see2 <- em_see[,c(1,3)]
      colnames(em_see2) <- c("master_list", "e0")
      write.csv(em_see2, "em_see.csv")
    } else {
      counter <- counter - 1
    }
    
    y <- y +1
  }

  while ( counter != 0 ) {
    
    if ( y > 10000 ) {
      y <- 1
      counter <- num_sim
    }
    
    print(paste("trial e", y, " counter ", counter, sep = ""))
    
    em_see <- read.csv("em_see.csv")[,-1]
    nprobes <- nrow(em_see)
    num_e0_probes <- sum(em_see$e0)
    
    nprobes <- nrow(em_see)
    num_e0_probes <- sum(em_see$e0)
    
    set.seed(floor(runif(1, min=0, max=(y*10))))
    neg_probes = floor(runif(1, min=0, max=(num_e0_probes/2)))
    pos_probes = floor(runif(1, min=0, max=(min((num_e0_probes/2),(nprobes-num_e0_probes)))))
    em_see$new <- em_see$e0
    
    # find indices of 0s and 1s
    zero_indices <- which(em_see$new == 0)
    one_indices <- which(em_see$new == 1)
    
    # sample and update for adding probes
    if (length(zero_indices) >= pos_probes && pos_probes > 0) {
      indices_to_update <- sample(zero_indices, pos_probes)
      em_see$new[indices_to_update] <- 1
    }
    
    # sample and update for removing probes
    if (length(one_indices) >= neg_probes && neg_probes > 0) {
      indices_to_update <- sample(one_indices, neg_probes)
      em_see$new[indices_to_update] <- 0
    }
    
    new <- em_see %>%
      select(master_list, new) %>%
      mutate(new = if_else(new == 0, NA_real_, new)) %>%
      na.omit() %>%
      select(-new) %>%
      setNames("Name")
    
    mx <- merge(x = new, y = miRNA, by.x = "Name", by.y = "...1") %>%
      select(-1)
    
    set.seed(42)
    clusters <- tibble(sample = colnames(mx), 
                       cluster = kmeans(t(mx), 3, iter.max = 10000)$cluster)
    
    ref <- read.csv("gdc_join_clinical.csv")
    notes <- merge(x = ref, y = clusters, by.x = "miRNA_filename", by.y = "sample") %>%
      select(-c(5,6,15:28))
    
    calculate_scores <- function(df, cluster_id) {
      df <- df %>%
        filter(cluster == cluster_id, Sample_Type == "Primary Tumor") %>%
        mutate(live_score = as.numeric(days_to_last_follow_up) / 624,
               death_score = as.numeric(days_to_death) / 761)
    }
    
    
    mxcls <- lapply(1:3, calculate_scores, df = notes)
    
    nx <- data.frame(samples = vapply(mxcls, nrow, numeric(1)),
                     dead = vapply(mxcls, function(df) sum(df$vital_status == "Dead"), numeric(1)),
                     live_weight = vapply(mxcls, function(df) sum(df$live_score, na.rm = TRUE), numeric(1)),
                     dead_weight = 0.01)
    
    total_dead <- sum(nx$dead)
    
    nx <- nx %>%
      rownames_to_column("cluster") %>%
      mutate(dead = if_else(dead == 0, 0.01, dead),
             dead_weight = if_else(dead > 0, if_else(vapply(mxcls, function(df) sum(df$death_score, na.rm = TRUE), numeric(1)) == 0, 0.01, vapply(mxcls, function(df) sum(df$death_score, na.rm = TRUE), numeric(1))), dead_weight),
             ratio = dead / samples,
             precision = dead / (dead + live_weight),
             recall = dead / total_dead,
             final_score = if_else((precision + recall) == 0, 0, 2 * (precision * recall) / (precision + recall)),
             rscore = ratio * dead_weight)
    
    nxa <- nx %>% arrange(final_score)
    nxb <- nx %>% arrange(rscore)
    
    na_best <- read.csv("na_best.csv") %>% select(-1)
    nb_best <- read.csv("nb_best.csv") %>% select(-1)
    
    if (round(nxa$final_score[3], 5) > round(na_best$final_score[3], 5) |
        (round(nxa$final_score[3], 5) >= round(na_best$final_score[3], 5) & round(nxb$rscore[1], 7) < round(nb_best$rscore[1], 7)) |
        (round(nxa$final_score[3], 5) >= round(na_best$final_score[3], 5) & round(nxb$rscore[1], 7) <= round(nb_best$rscore[1], 7) & nrow(new) > num_e0_probes)) {
      write.csv(nxa, "na_best.csv")
      write.csv(nxb, "nb_best.csv")
      
      print(paste("e", y, " is best", sep = ""))
      file.remove("em_see.csv")
      em_see2 <- em_see[,c(1,3)]
      colnames(em_see2) <- c("master_list", "e0")
      write.csv(em_see2, "em_see.csv")
    } else {
      counter <- counter - 1
    }
    
    y <- y +1
  }
  
#########################################################################  
  
  folder_name <- paste0("conv", i)
  if (!dir.exists(folder_name)) {
    dir.create(folder_name)
  }
  
  files_to_copy <- c("em_see.csv", "na_best.csv", "nb_best.csv", "trainingData.csv")
  
  folder_name <- paste0("conv", i)
  for (file_name in files_to_copy) {
    file.copy(file_name, file.path(folder_name, file_name))
  }
  cat("Files copied successfully!\n")
  
  # Save the current R environment
  save.image(file = "my_environment.RData")
  # load(file = "my_environment.RData")
  
  see <- read.csv("em_see.csv")
  data <- read_csv("trainingData.csv")
  filenames <- read.csv("filename_trainning.csv")
  
  see2 <- subset(see, e0 == 0)
  see2 <- see2[,-1]
  see2$e0 <- 1
  
  print("reduce hyperparameters")
  data2 <- merge(x=see2, y=data, by.x="master_list", by.y="...1", all = FALSE)
  rownames(data2) <- data2[,1]
  data2 <- data2[,-c(1,2)]
  colnames(data2) <- filenames$miRNA_filename
  
  write.csv(see2, "em_see.csv")
  write.csv(data2, "trainingData.csv")
  
}

