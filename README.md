#### The repository contains R codes and an example empirical multiple sequence <br /> alignmnet (MSA) for which the optimal model is selected uisng ModelTamer.

# Taming the selection of optimal substitution models in Phylogenomics

# ModelTamer 
ModelTamer is an automated tool for generating small subsamples of sites, expand the subsample to match the <br /> dimension of the full MSA, and find the optimal model of substitution (manuscript in review).
<br />

## Directory Structure 

The "code" directory contains files for two R function for: ``ModelTamer automatic`` (MT_automatic.R). This directory also contains executable IQTREE (iqtree) for analyzing MSA.
```
Runme(MT_automatic).R file has the source code for performing the ModelTamer analysis to select the best-fit model for the provided MSA.
```
The "data" directory contains one example empirical dataset which has been analyzed in the original article.  
```
Plants.fas          : This emprical dataset has 16 platns sequences and 4,246,454 sites with 190,352 distinct site patterns. The otimal model of substitution is GTR+G+R which is determined by analyzing the full MSA. 

```

## Introduction
The aoutomated ModelTamer tool:. 
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
All essential steps for ModelTamer analyses for the example empirical dataset can be performed in the  ``CodeOcean Capsule`` using the ``Reproducible Run`` functionality. 

## Reproducible Run
1. Set the ``Runme(MT_automatic).R`` as the file to run in CodeOcean. 
2. Click on the ``Reproducible Run`` button to perform the ModelTamer analysis.

## Github Repository
All the R codes and other instructions essential for ModelTamer analyses on a local computer are provided in the [Github Repository](https://github.com/ssharma2712/ModelTamer).
