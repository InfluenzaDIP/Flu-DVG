## Overview

This pipeline processes and analyzes single-cell RNA-seq (scRNA-seq) data using a series of bioinformatics tools, including the Cellranger toolkit, Seurat, Harmony, and others. The workflow is designed to ensure high-quality preprocessing, clustering, dimensionality reduction, and downstream analyses such as differential expression and functional enrichment. The provided code and instructions can reproduce the analysis demonstrated in the study.

---

## Data

The scRNA-seq data analyzed in this pipeline were deposited in the **NCBI Gene Expression Omnibus (GEO)** under the accession number **GSE277095**. These datasets can be used as a demo for reproducing the pipeline results.

---

## Requirements

### System
- **System architecture**: `x86_64, linux-gnu`

### Software Versions
- **R language**: R version `4.1.3`
- **Key tools**:  
  - Cellranger toolkit `v6.056`  
  - Seurat package `v4.057`  
  - Harmony algorithm `v0.158`  
  - clusterProfiler package `v3.1459`  
  - BLASTn `v2.1660`

---

## Installation

### Cellranger Toolkit

Please refer to the official Cellranger documentation for installation instructions:  
[Cellranger Installation Guide](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation)

### Seurat Package

Install the Seurat package in R using the following commands, or visit the official Seurat website:  
[Seurat Documentation](https://satijalab.org/seurat/)

```r
install.packages("Seurat")
```
