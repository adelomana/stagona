# Introduction
In this session we will apply unsupervised and supervised learning algorithms to single cell transcriptome data.
Please make sure you have your our required R packages installed.
Start RStudio and load the script _requisites.r_.

## Case 1. Unsupervised learning of melanoma malignant cells.
In our first example, we will study inter patient variability by learning cell clusters of malignant cells for different patients.
We will use data from a single-cell RNA sequencing (RNA-seq) study of 4,645 single cells isolated from 19 patients, profiling malignant, immune, stromal, and endothelial cells ([Tirosh et al., 2016](https://www.ncbi.nlm.nih.gov/pubmed/?term=27124452)).
We will use t-SNE (dimensionality reduction) and mean-shift clustering (unsupervised learning) to observe the differences between malignant cells from different patients.

## Case 2. Supervised learning of glioblastoma associated immune cell types.
In our second example, we will study intra patient variability by learning cell type of the non malignant cells for different patients.
We will use data from a single-cell RNA-seq study of We used single-cell RNA sequencing (RNA-seq) to profile 430 cells from five primary glioblastomas  ([Patel et al., 2014](https://www.ncbi.nlm.nih.gov/pubmed/24925914)).
We will use already known transcriptional signature to classify (supervised learning) immune cell types found in glioblastoma. 
