MT_automatic <- function(data_path, data_type = c("DNA", "AA"), Redo = FALSE, max.iter = 2){
  
  # data_path   : path for the sequence alignment in fasta format
  # data_type   : DNA or Amino Acid (AA)
  
  ########## Package required ##########
  
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  if (!requireNamespace("Biostrings", quietly = TRUE))
    BiocManager::install("Biostrings")
  
  if (!requireNamespace("stringr", quietly = TRUE))
    install.packages("stringr")
  
  ########## Library required ##########
  
  if (!library('Biostrings',logical.return = TRUE)){
    stop("'Biostrings' package not found, please install it to run MT_automtic")
  }
  
  if (!library('stringr',logical.return = TRUE)){
    stop("'stringr' package not found, please install it to run MT_automatic")
  }
  
  
  ######## ModelTamer Sampling #####################
  
  MT_sampler_auto <- function(MSA, num_of_dsp, distinct_positions, g, SU_name, sln = sln){
    
    # data_path: path for the sequence alignment in fasta format
    # g        : % of distinct site patterns required in a subsample
    # s        : Number of subsamples
    # r        : Number of upsamples
    
    
    directory <- getwd()
    SU_name <- SU_name           # subsample-upsample dataset name
    setwd(directory)
      
    expected_dsp <- ceiling(num_of_dsp*(g/100))
    initial_sample <- sample(distinct_positions, expected_dsp, replace = F, prob = NULL)
    expected_site <- ceiling((expected_dsp/length(unique(initial_sample)))*expected_dsp)
      
    sln1 <- sample(1:sln, expected_site, replace = F, prob = NULL)
    sln2 <- sample(sln1, sln, replace = T)
        
    print(c("% Distinct site pattern sampled = ", paste(signif((length(unique(sln1))/num_of_dsp)*100, 4),"%", sep ="") ))
    SU_dataset <- endoapply(motherfile, function(x) x[sln2])
        
    file_name <- paste(SU_name, sep = "")
    Biostrings::writeXStringSet(SU_dataset, file_name)
    
  }
      
  
  ##################################################
  
  f_name <- data_path                                                  # Mother MSA data path
  motherfile <- Biostrings::readAAStringSet(f_name, format = "fasta")  # Reading the mother file
  sln <- as.numeric(fasta.seqlengths(f_name)[1])                       # sequence length
  a <- as.data.frame(t(as.matrix(motherfile)))                         # convert MSA into a matrix
  
  a2 <- unique(a)                                                      # Finding distinct site patterns
  a2$ID <- as.vector(as.numeric(rownames(unique(a))))                  # Positions containing distinct site patterns
  a2 <- merge(a, a2)                                  
  
  distinct_positions <- a2$ID
  distinct_prob      <- prop.table(distinct_positions) 
  num_of_dsp <- length(unique(distinct_positions))                     # Number of distinct site patterns
  
  rm(a)
  rm(a2)  
  
  SU_name_base <- str_replace(basename(data_path), ".fasta", "")
  
  if(data_type == "DNA"){
    g_est <- floor((4594.2*num_of_dsp^(-1.043))*100)
  }else{
    g_est <- floor(((4/20)*4594.2*num_of_dsp^(-1.043))*100)
  }
 
  
  if(g_est >= 63){
    
    print("ModelTamer will analyze 63% of distinct site patterns")
    print("ModelTamer suggests to analyze full MSA")
    
    g_est <- 100
  }else{
    g_est <- g_est
  }

  if(Redo == "FALSE"){
    
    SU_name1 <- paste(SU_name_base, "_g_", g_est, ".fasta", sep = "")
    MT_sampler_auto(motherfile, num_of_dsp = num_of_dsp, distinct_positions = distinct_positions,
                    g = g_est, SU_name = SU_name1, sln = sln)
    
    if(as.character(Sys.info()[1]) == "Windows"){
      ex_cmd <- paste("iqtree -s", SU_name1, "-m MF -nt 1 -quiet",  sep = " " )
      system(ex_cmd)
    }else{
      ex_cmd <- paste("./iqtree -s", SU_name1, "-m MF -nt 1 -quiet",  sep = " " )
      system(ex_cmd)
    }
    
   
    sink("MT_output.txt", append = T)
    
    a_log <- readLines(paste(SU_name1, ".log", sep = ""))
    cline1 <- which(a_log == grep(pattern = "Best-fit model:", a_log, value = TRUE))
    a_log <- a_log[cline1]
    final_model <- sapply(a_log, function(x){unlist(strsplit(x, "\\s+"))})[,][3]
    print(paste("Final Model (BIC):", final_model, sep = " "))
    sink()
    
  }else{
    g1 <- (1:max.iter)*g_est
    g1 <- c(g1[g1<63], 100)
    g1 <- g1[1:max.iter]
    
    for(i in 1:length(g1)){
      SU_name1 <- paste(SU_name_base, "_g_", g1[i], ".fasta", sep = "")
      MT_sampler_auto(motherfile, num_of_dsp = num_of_dsp, distinct_positions = distinct_positions,
                      g = g1[i], SU_name = SU_name1, sln = sln)
      
      if(as.character(Sys.info()[1]) == "Windows"){
        ex_cmd <- paste("iqtree -s", SU_name1, "-m MF -nt 1 -quiet",  sep = " " )
        system(ex_cmd)
      }else{
        ex_cmd <- paste("./iqtree -s", SU_name1, "-m MF -nt 1 -quiet",  sep = " " )
        system(ex_cmd)
      }  # else
      
      sink("MT_output.txt", append = T)
      print(c("Iteration:", i))
      cat("\n\n")
      a_log <- readLines(paste(SU_name1, ".log", sep = ""))
      cline1 <- which (a_log == grep(pattern = "Best-fit model:", a_log, value = TRUE))
      a_log <- a_log[cline1]
      final_model <- sapply(a_log, function(x){unlist(strsplit(x, "\\s+"))})[,][3]
      print(paste("Final Model (BIC) in iteration", i, ":", final_model, sep = " "))
      cat("\n\n\n")
      sink()
      
    }    # for loop
  }      # else

  

}


###############  MT_aggregator  #######################

MT_aggregator <- function(log_file_path, iqtree_file_path, data_t = c("DNA", "AA")){
  
  # log_file_path           : path for the log file 
  # iqtree_file_path        : % of distinct site patterns required in a subsample
  # data_type               : DNA or Amino Acid (AA)
  
  #############log to Substitution matrix extraction ####################
  log2submat <- function(iqtree_file){
    cline1 <- which (iqtree_file == grep(pattern = "Rate matrix Q:", iqtree_file, value = TRUE))
    iqtree_file <- iqtree_file[(cline1+2):(cline1+5)]
    iqtree_file <- sapply(iqtree_file, function(x){unlist(strsplit(x, "\\s+"))})
    
    iqtree_file <- data.frame(matrix(unlist(iqtree_file), nrow=4, byrow=TRUE),stringsAsFactors=FALSE)
    rownames(iqtree_file) <- iqtree_file[,2]
    iqtree_file <- iqtree_file[,-c(1,2)]
    colnames(iqtree_file) <- rownames(iqtree_file)
    return(iqtree_file)
  }
  ##########################################################################
  
  output <- list()
  
  a_iqtree <- readLines(iqtree_file_path)
  a_log    <- readLines(log_file_path)
  
  
  cline2 <- which (a_log == grep(pattern = "Corrected Akaike Information Criterion:", a_log, value = TRUE))
  
  
  output[[1]] <- a_log[cline2+2]
  
  output[[2]] <- a_log[cline2-1]
  output[[3]] <- a_log[cline2]
  output[[4]] <- a_log[cline2+1]
  
  cline1 <- which (a_iqtree == grep(pattern = "Input data:", a_iqtree, value = TRUE))
  output[[5]] <- a_iqtree[cline1]
  output[[6]] <- a_iqtree[cline1+1]
  output[[7]] <- a_iqtree[cline1+2]
  output[[8]] <- a_iqtree[cline1+3]
  output[[9]] <- a_iqtree[cline1+4]
  
  
  cline2 <- which (a_log == grep(pattern = "NOTE: ModelFinder requires", a_log, value = TRUE))
  output[[10]] <- a_log[cline2]
  
  cline2 <- which (a_log == grep(pattern = "CPU time for ModelFinder:", a_log, value = TRUE))
  output[[11]] <- a_log[cline2]
  output[[12]] <- a_log[cline2+1]
  
  sink("MT_output.txt")
  print(unlist(output))
  sink()
  
  if(data_t == "DNA"){
    
    submat <- log2submat(a_iqtree)
    submat <- sapply(submat, as.numeric)
    rownames(submat) <- colnames(submat)
    output[[13]] <- "Substitution Matrix (Q)"
    sink("MT_output.txt", append = T)
    cat("\n\n")
    print(unlist(output[[13]]))
    print(submat)
    sink()
  }else{
    output <- output
  }
  
}

