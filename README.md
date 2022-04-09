#### The repository contains R codes and an example empirical multiple sequence alignmnet (MSA) for which the optimal model is selected uisng ModelTamer.

# Taming the selection of optimal substitution models in Phylogenomics

# ModelTamer 
``ModelTamer`` is an automated tool for generating small subsamples of sites, expand the subsample to match the dimension of the full MSA, and find the optimal model of substitution (manuscript in review).
<br />

## Directory Structure 

The "Code" directory contains files for three R functions:
```
MT_sampler.R    : Create subsample and upsample dataset (SU dataset) using user specified percentage of distinct site patterns.
MT_aggregator.R : Extract analysis information from IQTREE output. 
MT_automatic.R  : An aumtomated tool for performing ModelTamer. 
```
The "Example" directory contains one example empirical dataset which has been analyzed in the original article.  
```
Plants.fas      : This emprical dataset has 16 platns sequences and 4,246,454 sites with 190,352 distinct site patterns. The otimal model of substitution is GTR+G+R which is determined by analyzing the full MSA. 

```

## Introduction
The aoutomated ModelTamer tool:. 
<br />

```
MT_sampler(data_path, g, s=1, r =1)

data_path         : input sequence alignment in fasta format that will be used for little bootstrap analyses. 
g                 : a numeric input for percentage of distinct site patterns we want to sample in the SU dataset.
s                 : a numeric input for the number of subsamples we want to analyze. 
r                 : number of upsample dataset per subsample. It is equal to 1 (r = 1) by default. 
```
<br />

```
MT_aggregator(log_file_path, iqtree_file_path, data_t = c("DNA", "AA"))

data_path         : input sequence alignment in fasta format that will be used for little bootstrap analyses. 
g                 : percentage of distinct site patterns we want to sample in the SU dataset.
data_t            : a character input argument to specify the base code type whether the MSA contains DNA ("DNA") or amino acid ("AA") sequences. 

```
<br />

```
MT_automatic(data_path, data_type = c("DNA", "AA"), Redo = FALSE, max.iter = 2)

data_path         : input sequence alignment in fasta format that will be used for little bootstrap analyses. 
data_type         : a character input argument to specify the base code type whether the MSA contains DNA ("DNA") or amino acid ("AA") sequences.
Redo              : a logical input argument to specify the program to re-check the estimated otimal model found in first step by increasing the subsample size. 
max.iter          : a numeric input argument that specifies the number of steps to be performed by ModelTamer by increasing the subsample size.

```

ModelTamer uses ModelFinder from IQTREE by default for performing calculating the Maximum Likelihood (ML) fit for each of the model tested. The IQTREE software can be downloaded from http://www.iqtree.org/. The ModelTamer is designed to work on both both Linux and Windows operating system. However, it is recomened to used the operating system specific version of IQTREE for performing ModelTamer. 

<br />

## Getting Started with ModelTamer

To perform the model selection analysis using ModelTamer on your local computer, please follow these steps:<br /><br />
1.	Download and install R (https://www.r-project.org/) and Rstudio (https://rstudio.com/products/rstudio/download/).<br />
2.	Download the github repository containing R codes and the example dataset on the local computer. <br />
3.	In the Rstudio session, type ``setwd(“directory path”)`` to change the working directory to the folder that contains ``MT_sampler.R``, ``MT_aggregator.R``, and and ``MT_automatic.R`` function<br />
4.	Type ``source(lb_sampler.R)``, and ``source(lb_aggregator.R)`` or  ``source(lb_precision.R)`` to make these functions available in the global environment. <br />
5.	Download and install an ML tree inference software (e.g., IQ-TREE). <br />
6.	Install the following R packages if thay are not installed. 

```R
install.packages("BiocManager")
BiocManager::install("Biostrings")
install.packages("stringr")
``` 

<br />

## ModelTamer Analyses for the Example Dataset:

<br />
To create SU dataset for the example Plants MSA, please follow these steps:<br /><br />
1.	Download the ``Example`` directory on the local computer. <br />
2.	Run the function in the R session by typing 

```R
MT_sampler("~/Example/Plants.fas", g = 0.1, s = 1, r = 1)
```

This function will create a directory in the working directory:

```
MT_Subsample_g_0.1
```

The directory will contain one SU dataset. For example 

```
Plants_sub_1up_1.fasta

```
3.	Find the best-fit model for the given SU dataset by using:<br />

``` 
iqtree -s ~/Example/MT_Subsample_g_0.1/Plants_sub_1up_1.fasta -m MF
```

The .log and .iqtree file will be restored in the ``MT_Subsample_g_0.1`` folder like Plants_sub_1up_1.fasta.log and Plants_sub_1up_1.fasta.iqtree. <br /><br />
4.	For the final step, type 

```R
MT_aggregator("~/Example/MT_Subsample_g_0.1/Plants_sub_1up_1.fasta.log", "~/Example/MT_Subsample_g_0.1/Plants_sub_1up_1.fasta.iqtree", data_t = "DNA" )
```
The function will output a text file with ModelTamer analysis information, and the name of the output will be `` MT_output.txt``.<br />

<br />

## Automated ModelTamer Analyses for the Example Dataset:

```R
MT_automatic("~/Example/Plants.fas", data_type = "DNA", Redo = FALSE, max.iter = 2 )
```
This function will select model for the Plants dataset automatically and output ``MT_output.txt`` which contains the final model selected by ``ModelTamer``.

#### Software and Packages' Version:

<br />

All R codes were tested using R version 3.6.3 in R studio (version 1.2.5033).
<br />  
R packages used:
<br />

```
-BiocManager (version 1.30.10)
-Biostrings  (version 2.54.0)
-stringr     (version 1.4.0)
```

<br />
<br />
