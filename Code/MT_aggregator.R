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
  
  cline1 <- which (a_iqtree == grep(pattern = "Input data:", a_iqtree, value = TRUE))
  cline2 <- which (a_log == grep(pattern = "Corrected Akaike Information Criterion:", a_log, value = TRUE))
  
  output[[1]] <- "\n"
  output[[2]] <- "Best-fit model by Information Criteria (ModelTamer)"
  output[[3]] <- "\n\n"
  output[[4]] <- a_log[cline2+1]
  output[[5]] <- "\n"
  output[[6]] <- a_log[cline2-1]
  output[[7]] <- "\n"
  output[[8]] <- a_log[cline2]
  
  output[[9]] <- "\n\n"
  output[[10]] <- a_iqtree[cline1]
  output[[11]] <- "\n"
  output[[12]] <- paste(a_iqtree[cline1+3], "(reported by IQTREE)", sep = " ")
  output[[11]] <- "\n"
  output[[12]] <- paste("Distinct patterns used by ModelTamer: ",sapply(a_iqtree[cline1+4], function(x){unlist(strsplit(x, "\\s+"))})[6], sep = "")
  output[[13]] <- "\n\n"
  
  cline2 <- which (a_log == grep(pattern = "NOTE: ModelFinder requires", a_log, value = TRUE))
  a <- sapply(a_log[cline2], function(x){unlist(strsplit(x, "\\s+"))})[4:6]
  output[[14]] <- paste("Peak memory used by ModelTamer: ",a[1], a[2], a[3], sep = " ")
  output[[15]] <- "\n"
  
  cline2 <- which (a_log == grep(pattern = "CPU time for ModelFinder:", a_log, value = TRUE))
  a <- sapply(a_log[cline2], function(x){unlist(strsplit(x, "\\s+"))})[5:7]
  output[[16]] <- paste("CPU time for ModelTamer: ",a[1], a[2], a[3], sep = " ")
  
  output[[17]] <- "\n"
  a <- sapply(a_log[cline2+1], function(x){unlist(strsplit(x, "\\s+"))})[5:7]
  output[[18]] <- paste("Wall-clock time for ModelTamer: ",a[1], a[2], a[3], sep = " ")
  
  
  sink("MT_output.txt")
  cat(unlist(output))
  sink()
  
  if(data_t == "DNA"){
    
    submat <- log2submat(a_iqtree)
    submat <- sapply(submat, as.numeric)
    rownames(submat) <- colnames(submat)
    
    sink("MT_output.txt", append = T)
    cat("\n\n")
    cat("Substitution Matrix (Q)")
    cat("\n")
    print(submat)
    sink()
  }else{
    output <- output
  }
  
}
