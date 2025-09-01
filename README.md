# Spatial Transcriptomics Analysis – Visium HD Mouse Pancreas TMAs

This repository contains a comprehensive pipeline for analyzing spatial transcriptomics data from four Visium HD tissue microarrays (TMAs) of mouse pancreatic islets. The study compares PolgMUT vs wild-type mice, focusing on differential gene expression and pathway changes in endocrine cell populations.

## Overview

The pipeline processes 4 TMAs with pancreatic islet data, performs batch correction using Harmony, identifies pure endocrine cell types, and conducts differential expression analysis between genotypes.

**Key Features:**
- Automated SpaceRanger processing for all TMAs
- Coordinate matching between old and new data formats
- Harmony batch correction across TMAs
- Pure cell type identification using hormone markers
- Comprehensive differential expression analysis
- Publication-ready visualizations

---

## Prerequisites

### Software Requirements
- **SpaceRanger**: 2.1.1 (automatically downloaded by setup script)
- **R**: ≥4.2.0 with required packages (automatically installed)
- **System**: Linux/macOS with ≥30GB RAM, ≥14 CPU cores recommended

### R Package Dependencies
The pipeline automatically installs these packages if missing:
- Seurat (≥4.0)
- harmony
- ggplot2, patchwork, dplyr, tidyr
- RColorBrewer, Matrix, arrow, hdf5r

---

## Quick Start

### 1. Setup SpaceRanger and Reference Genome
```bash
# Download and setup SpaceRanger + mouse reference genome
chmod +x setup_spaceranger.sh
./setup_spaceranger.sh
```

### 2. Configure and Run SpaceRanger
```bash
# Edit paths in the script to match your data location
nano run_spaceranger_all.sh

# Update these paths:
# - TRANSCRIPTOME path
# - Library CSV files
# - H&E image paths  
# - Loupe alignment files

# Run SpaceRanger for all TMAs
chmod +x run_spaceranger_all.sh
./run_spaceranger_all.sh
```

### 3. Run Spatial Analysis Pipeline
```R
# Load the analysis functions
source("polg_spatial_analysis.R")

# Run complete pipeline
results <- run_polg_analysis_pipeline(
  hd_data_dir = "/path/to/spaceranger/outputs",
  old_data_dir = "/path/to/islet/annotations", 
  output_dir = "polg_analysis_results"
)
```

---

## Data Structure

### Required Input Files

```
project/
├── spaceranger_outputs/           # SpaceRanger HD outputs
│   ├── TMA1/
│   │   └── binned_outputs/
│   │       └── square_008um/      # 8μm resolution data
│   ├── TMA2/
│   ├── TMA3/
│   └── TMA4/
└── islet_annotations/             # Legacy islet annotations
    ├── TMA1/
    │   ├── SpatialCoords.csv      # Original spatial coordinates
    │   ├── Islet_POLG_TMA1.csv    # POLG islet annotations
    │   └── Islet_WT_TMA1.csv      # WT islet annotations
    ├── TMA2/
    ├── TMA3/
    └── TMA4/
```

### SpaceRanger Configuration Files

Each TMA requires a `libraries.csv` file:
```csv
fastqs,sample,library_type
/path/to/fastq/TMA1,TMA1_sample,Gene Expression
```

---

## Analysis Pipeline Details

### Step 1: Data Loading and Coordinate Matching
- Loads Visium HD data at 8μm resolution
- Matches coordinates between new HD format and legacy annotations
- Creates islet-specific Seurat objects with genotype labels

### Step 2: Batch Correction with Harmony
- SCTransform normalization across all TMAs
- PCA dimensionality reduction
- Harmony batch correction to integrate TMAs
- UMAP visualization of corrected data

### Step 3: Cell Type Identification
Cell type identification incorporates comprehensive gene-specific marker panels and hormone expression thresholds. The pipeline uses established pancreatic cell markers to classify pure cell populations and identify tissue architecture:

**Primary hormone markers for pure endocrine cells:**
- **Ins1, Ins2** (Insulin 1 & 2): β-cell markers
- **Gcg** (Glucagon): α-cell marker  
- **Sst** (Somatostatin): δ-cell marker
- **Ppy** (Pancreatic polypeptide): γ-cell marker

**Additional cell type-specific markers:**
- **β-cells**: Pdx1, MafA, Nkx6.1, Slc2a2 (Glut2), Ucn3
- **α-cells**: Arx, Irx1, Irx2, Mafb, Tm4sf4
- **δ-cells**: Hhex, Pdx1 (low expression), Lepr
- **γ-cells**: Mafb, Ghrl (ghrelin)
- **Ductal cells**: Krt19, Sox9, Hnf1b, Cftr, Spp1
- **Acinar cells**: Amy1, Amy2a, Cpa1, Ptf1a, Ctrb1
- **Endothelial cells**: Pecam1, Cdh5, Tie1, Plvap
- **Immune cells**: Ptprc (Cd45), Cd68, Cd3e, Lyz2

**Classification criteria for pure endocrine cells:**
- **β-cells**: High Ins1/Ins2 (>0.5), low Gcg/Sst (<0.5)
- **α-cells**: High Gcg (>0.5), low Ins1/Ins2/Sst (<0.2-0.5)  
- **δ-cells**: High Sst (>0.5), low Ins1/Ins2/Gcg (<0.2-0.5)
- **Mixed Endocrine**: Co-expression of multiple hormones
- **Other**: Non-endocrine or unidentified cells

Only "pure" cell types with clear single-hormone dominance are used for downstream differential expression analysis, ensuring cell-type-specific comparisons between POLG and WT genotypes.

### Step 4: Differential Expression Analysis
- Wilcoxon rank-sum tests comparing POLG vs WT
- Separate analysis for each pure cell type
- Multiple testing correction (Bonferroni)
- Minimum 10 cells per group per comparison

### Step 5: Visualization and Export
- Volcano plots for each cell type
- UMAP plots colored by genotype and cell type
- Batch correction comparison plots
- All results exported as CSV and PNG files


## Key Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| Resolution | 8μm | Spatial resolution for analysis |
| Min cells per group | 10 | Minimum cells for DE analysis |
| Coordinate tolerance | 10 pixels | Matching tolerance for coordinates |
| Harmony dims | 1:30 | PCA dimensions for batch correction |
| DE test | wilcox | Statistical test for DE analysis |
| FC threshold | 0.1 | Log2 fold-change threshold |
| P-value threshold | 0.05 | Adjusted p-value significance |


---

## Contact

**Alana Mullins**  
Email: a.mullins2@newcastle.ac.uk





