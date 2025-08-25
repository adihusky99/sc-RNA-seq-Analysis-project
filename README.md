# Overview
This project performs comprehensive exploratory data analysis (EDA) on RNA-sequencing data from central nervous system (CNS) cancer cell lines from the Cancer Cell Line Encyclopedia (CCLE). The analysis identifies distinct molecular subtypes and characterizes their biological signatures through unsupervised machine learning approaches.

## Dataset

Source: Cancer Cell Line Encyclopedia (CCLE)
Samples: 63 CNS cancer cell lines
Features: ~19,000 protein-coding genes
Data Type: Raw RNA-seq count data with associated metadata

## Files Structure

data/
├── CCLE_RNAseq_genes_counts_20180929.gct.gz    # Raw count matrix
├── Cell_lines_annotations_20181226.txt          # Sample metadata
outs/
├── project_ccle_counts.csv                      # Processed count matrix
├── project_ccle_meta.csv                        # Processed metadata
├── project_ccle_counts_subset.csv               # CNS-specific counts
└── project_ccle_meta_subset.csv                 # CNS-specific metadata


## Key Analyses
### 1. Data Preprocessing

Gene filtering: Retained only protein-coding genes using biomaRt
Duplicate handling: Merged duplicated gene entries by summing counts
Quality control: Filtered low-count genes (≥100 total counts)
Normalization: Applied variance stabilizing transformation (VST)

#### 2. Principal Component Analysis (PCA)

Feature selection: Top 500 most variable genes
Variance explained: PC1 (15%), PC2 (10%)
Key findings:

PC1: Driven by collagen genes (COL1A1, COL1A2, COL3A1, COL5A1, COL6A2, COL6A3)
PC2: Characterized by neural/glial markers (GFAP, SOX2, GPM6B, GPC4, CNR1, COL11A1)



#### 3. Hierarchical Clustering

Method: Complete linkage with Euclidean distance
Optimal clusters: k=3 (determined through systematic evaluation)
Cluster sizes: 37, 19, and 7 samples

#### 4. Differential Expression Analysis

Method: DESeq2 with negative binomial modeling
Significance threshold: padj < 0.001
Top differentially expressed genes identified between clusters

### Main Results
Three Distinct Molecular Subtypes Identified:

#### Cluster 1 (n=37): ECM/Stromal-like

High expression of collagen genes
Enriched for extracellular matrix remodeling pathways
Represents fibroblast-like or stromal characteristics


#### Cluster 2 (n=19): Neural/Glial-like

High expression of neural differentiation markers
Enriched for CNS-specific transcriptional programs
Represents more differentiated neural phenotypes


#### Cluster 3 (n=7): Intermediate

Balanced expression of both gene sets
Mixed or transitional molecular state



## Dependencies
R Packages Required:
r# Core analysis
library(DESeq2)
library(biomaRt)

##### Data manipulation
library(dplyr)
library(tidyr)

##### Visualization
library(ggplot2)
library(ComplexHeatmap)
library(EnhancedVolcano)
library(patchwork)
library(scales)
library(ggrepel)
Installation:
r# Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("DESeq2", "biomaRt", "ComplexHeatmap", "EnhancedVolcano"))

##### CRAN packages
install.packages(c("dplyr", "tidyr", "ggplot2", "patchwork", "scales", "ggrepel"))


## Key Visualizations Generated

#### Quality Control Plots:

Mean vs. variance scatter plots (raw vs. VST)
Expression distribution comparisons


#### PCA Plots:

PC1 vs PC2 colored by various metadata
PC3 vs PC4 analysis
Elbow plot for variance explanation


#### Clustering Visualizations:

Hierarchical clustering heatmap with annotations
PCA overlay with cluster assignments


#### Differential Expression:

Volcano plots highlighting significant genes
Expression bubble plots in PC space
Boxplots of top differentially expressed genes



