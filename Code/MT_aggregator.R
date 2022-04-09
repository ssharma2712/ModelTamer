MT_aggregator <- function(log_file_path, iqtree_file_path, data_t = c("DNA", "AA")){
  
  # log_file_path           : path for the log file 
  # iqtree_file_path        : % of distinct site patterns required in a subsample
  # data_t                  : DNA or Amino Acid (AA)
  
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




