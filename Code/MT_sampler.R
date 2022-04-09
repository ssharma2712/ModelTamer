
MT_sampler <- function(data_path, g, s = 1, r = 1){
  
  # data_path: path for the sequence alignment in fasta format
  # g        : % of distinct site patterns required in a subsample
  # s        : Number of subsamples
  # r        : Number of upsamples
  
  
  ########## Package required ##########
  
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  if (!requireNamespace("Biostrings", quietly = TRUE))
    BiocManager::install("Biostrings")
  
  if (!requireNamespace("stringr", quietly = TRUE))
    install.packages("stringr")
  
  ########## Library required ##########
  
  if (!library('Biostrings',logical.return = TRUE)){
    stop("'Biostrings' package not found, please install it to run lb_sampler")
  }
  
  if (!library('stringr',logical.return = TRUE)){
    stop("'stringr' package not found, please install it to run lb_sampler")
  }
  ######################################
  
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
  
  directory <- getwd()
  SU_name <- str_replace(basename(data_path), ".fasta", "")           # subsample-upsample dataset name
  
  for(i in 1:s){
    setwd(directory)
    
    expected_dsp <- ceiling(num_of_dsp*(g/100))
    initial_sample <- sample(distinct_positions, expected_dsp, replace = F, prob = NULL)
    expected_site <- ceiling((expected_dsp/length(unique(initial_sample)))*expected_dsp)
    
    subup_folder <- paste("MT_Subsample_g_", g, sep = "") 
    dir.create(subup_folder)
    setwd(subup_folder)
   # sub_name <- paste(SU_name, i, "_sub_", sep = "")
    
    for(j in 1:r){
      
      sln1 <- sample(1:sln, expected_site, replace = F, prob = NULL)
      sln2 <- sample(sln1, sln, replace = T)
      
      print(c("% Distinct site pattern sampled = ", paste(signif((length(unique(sln1))/num_of_dsp)*100, 4),"%", sep ="") ))
      
      SU_dataset <- endoapply(motherfile, function(x) x[sln2])
      
      file_name <- paste(SU_name, "_sub_", i, "up_", j, ".fasta", sep = "")
      Biostrings::writeXStringSet(SU_dataset, file_name)
      print(c("Subsample =", i, "Upsample =", j))
      rm(SU_dataset)
    }
   
  }
  setwd(directory)
}


