# Spatial Transcriptomics Analysis – Visium HD Mouse Pancreas TMAs

This repository contains code and documentation for analysing spatial transcriptomics data from four Visium HD tissue microarrays (TMAs) of fixed mouse pancreas. The study compares wild-type and PolgMUT mice across both sexes, focusing on transcriptional and pathway changes in endocrine regions.

## Overview

Each TMA contains 12 tissue cores (total 48 cores). Three regions of interest (ROIs) per sample were identified using Loupe Browser to isolate endocrine areas. Data are processed using SpaceRanger and analysed in R with Seurat and pathway enrichment tools.

---

## 1. Data Preprocessing

### a. Run SpaceRanger

Use the `run_spaceranger_all.sh` script to run `spaceranger count` for each tissue core.

- **Inputs**: FASTQ files and matched H&E images  
- **Outputs**: Spatial gene expression matrices per core

### b. Generate Metadata

Create a metadata table including:

- Genotype (wild-type or PolgMUT)
- Sex (male or female)
- TMA ID
- Sample condition
- Region of interest (ROI) annotation

---

## 2. ROI Selection and Annotation

- Use Loupe Browser to annotate endocrine ROIs in each core  
- Export ROI barcodes and incorporate into the metadata  
- Subset the dataset using these endocrine-specific barcodes

---

## 3. Data Loading and Integration in R

### a. Load and Merge Data

- Load cores using `Load10X_Spatial()` in Seurat  
- Merge into one Seurat object  
- Add metadata (condition, sex, TMA ID, ROI)

### b. Normalisation and Dimensionality Reduction

- Normalise with SCTransform  
- Run PCA and UMAP for visualisation  
- Cluster the data using `FindClusters()`

---

## 4. Cell Type Annotation

- Assign cell types using known markers (e.g. beta, alpha, delta, ductal)  
- Label clusters and visualise spatial distribution across TMAs

---

## 5. Differential Expression and Cell Composition

### a. Differential Expression

- Compare wild-type vs Polg^mut within each cell type and sex  
- Use `FindMarkers()` to identify differentially expressed genes

### b. Cell Composition Analysis

- Quantify proportion of each cell type in ROIs per core  
- Compare cell type frequencies across conditions

---

## 6. Pathway Enrichment

Perform enrichment analysis on differentially expressed genes using:

- KEGG  
- MitoCarta  

---

## 7. Output and Visualisation

- UMAPs coloured by cell type, condition, and genotype  
- Heatmaps of marker and differentially expressed genes  
- Per-core plots showing endocrine cell composition  
- Export UMAP coordinates and cluster IDs for downstream analysis

---

## Requirements

- SpaceRanger: 2.1.0  
- R: 4.2.0 with Seurat, tidyverse, and relevant enrichment packages  
- Loupe Browser (for ROI annotation)  

---

## Contact

For questions or feedback, please contact:  
Alana Mullins




