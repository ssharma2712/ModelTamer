#### The repository contains R codes and an example empirical multiple sequence alignmnet (MSA) for which the optimal model is selected uisng ModelTamer.

# Taming the selection of optimal substitution models in Phylogenomics

# ModelTamer 
ModelTamer is an automated tool for generating small subsamples of sites, expand the subsample to match the dimension of the full MSA, and find the optimal model of substitution (manuscript in review).
<br />

## Directory Structure 

The "Code" directory contains files for three R functions:
```
MT_sampler.R    : Create subsample and upsample dataset (SU dataset) using user specified percentage of distinct site patterns.
MT_aggregator.R : Extract analysis information from IQTREE output. 
MT_automatic.R  : An aumtomated tool for performing ModelTamer. 
```
The "data" directory contains one example empirical dataset which has been analyzed in the original article.  
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

## Getting Started with the CodeOcean Capsule
<br />

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
install.packages("ape")
install.packages("phangorn")
``` 

## Reproducible Run
1. Set the ``Runme(MT_automatic).R`` as the file to run in CodeOcean. 
2. Click on the ``Reproducible Run`` button to perform the ModelTamer analysis.

## Github Repository
All the R codes and other instructions essential for ModelTamer analyses on a local computer are provided in the [Github Repository](https://github.com/ssharma2712/ModelTamer).
