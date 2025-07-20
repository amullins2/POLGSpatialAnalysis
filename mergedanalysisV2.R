# spatial transcriptomics analysis pipeline for mouse pancreatic tmas with comprehensive QC
#
# this script performs the following tasks:
# - installs and loads required cran and bioconductor packages
# - loads spatial transcriptomics data from multiple tmas with loupe cluster annotations
# - performs comprehensive quality control and filtering
# - merges and normalises data across tmas
# - visualises qc metrics and marker gene expression for endocrine cell types
# - performs sketched clustering and projects clusters onto full data
# - identifies and highlights endocrine cell clusters spatially
# - subsets endocrine cells and performs differential gene expression analysis comparing genotypes (polg_mut vs wt)
# - converts gene symbols to entrez ids for enrichment analyses
# - performs mitocarta and kegg pathway enrichment on degs
# - runs gsea using kegg gene sets from msigdb
# - saves processed combined object and deg results for downstream analyses
#
# author: alana mullins
# date: 18/05/2025
# =============================================================================

# install cran packages, if not there
cran_packages <- c("Seurat", "dplyr", "ggplot2", "patchwork", "gridExtra")
new_cran <- cran_packages[!(cran_packages %in% installed.packages()[,"Package"])]
if(length(new_cran)) install.packages(new_cran)

# install bioconductor packages, if not there
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

bioc_packages <- c("clusterProfiler", "org.Mm.eg.db", "msigdbr", "fgsea")
new_bioc <- bioc_packages[!(bioc_packages %in% installed.packages()[,"Package"])]
if(length(new_bioc)) BiocManager::install(new_bioc)

# load packages
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(gridExtra)
library(clusterProfiler)
library(org.Mm.eg.db)   
library(msigdbr)        
library(fgsea)          

# load metadata (fixed column names to match your CSV)
metadata <- read.csv("TMAMetadata.csv", stringsAsFactors = FALSE)
# Rename columns to match script expectations
colnames(metadata) <- c("sample_id", "core_id", "tma_id", "genotype", "sex")

print("Metadata structure:")
print(head(metadata))
print(table(metadata$genotype))

tma_dirs <- c(
  "/TMA1/",
  "/TMA2/",
  "/TMA3/",
  "/TMA4/"
)

loupe_csvs <- c(
  "/TMA1_loupe.csv",
  "/TMA2_loupe.csv",
  "/TMA3_loupe.csv",
  "/TMA4_loupe.csv"
)

tma_names <- c("TMA1", "TMA2", "TMA3", "TMA4")

# function to load and annotate one tma 
load_and_annotate_TMA <- function(tma_dir, loupe_csv_path, tma_name, metadata) {
  cat("Loading TMA:", tma_name, "\n")
  object <- Load10X_Spatial(data.dir = tma_dir, bin.size = c(8, 16))
  
  # loupe clusters (endocrine annotations)
  if(file.exists(loupe_csv_path)) {
    loupe_annot <- read.csv(loupe_csv_path)
    rownames(loupe_annot) <- loupe_annot$barcode
    common_barcodes <- intersect(Cells(object), rownames(loupe_annot))
    object$Loupe_Cluster <- NA
    object$Loupe_Cluster[common_barcodes] <- loupe_annot[common_barcodes, "cluster"]
    cat("Loupe annotations loaded for", length(common_barcodes), "spots\n")
  } else {
    warning("Loupe CSV not found:", loupe_csv_path)
    object$Loupe_Cluster <- NA
  }
  
  # extract core_id & add metadata
  object$core_id <- sapply(strsplit(Cells(object), "-"), `[`, 1)
  
  # Add metadata fields
  for (field in c("sample_id", "genotype", "sex")) {
    object[[field]] <- metadata[match(object$core_id, metadata$core_id), field]
  }
  
  object$TMA <- tma_name
  return(object)
}

# load all tmas
tma_list <- lapply(seq_along(tma_dirs), function(i) {
  load_and_annotate_TMA(tma_dirs[i], loupe_csvs[i], tma_names[i], metadata)
})
names(tma_list) <- tma_names

# =============================================================================
# STAGE 1: INDIVIDUAL TMA QC BEFORE MERGING
# =============================================================================

cat("\n=== STAGE 1: INDIVIDUAL TMA QC ===\n")

# Function to calculate QC metrics for each TMA
calculate_qc_metrics <- function(object, assay_name = "Spatial.008um") {
  DefaultAssay(object) <- assay_name
  
  # Calculate mitochondrial gene percentage
  object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^mt-")
  
  # Calculate ribosomal gene percentage
  object[["percent.ribo"]] <- PercentageFeatureSet(object, pattern = "^Rp[sl]")
  
  # Calculate hemoglobin gene percentage (tissue artifact)
  object[["percent.hb"]] <- PercentageFeatureSet(object, pattern = "^Hb[^(p)]")
  
  return(object)
}

# Apply QC metrics to all TMAs
tma_list <- lapply(tma_list, calculate_qc_metrics)

# QC visualization for each TMA
for(i in seq_along(tma_list)) {
  cat("QC for", names(tma_list)[i], "\n")
  
  # Print basic stats
  cat("Number of spots:", ncol(tma_list[[i]]), "\n")
  cat("Number of genes:", nrow(tma_list[[i]]), "\n")
  
  # QC violin plots
  DefaultAssay(tma_list[[i]]) <- "Spatial.008um"
  p1 <- VlnPlot(tma_list[[i]], features = c("nFeature_Spatial.008um", "nCount_Spatial.008um", 
                                           "percent.mt", "percent.ribo"), ncol = 2, pt.size = 0)
  print(p1)
  
  # Spatial QC plots
  p2 <- SpatialFeaturePlot(tma_list[[i]], features = c("nCount_Spatial.008um", "nFeature_Spatial.008um"), ncol = 2)
  print(p2)
}

# =============================================================================
# STAGE 2: MERGE TMAS AND COMPREHENSIVE QC
# =============================================================================

cat("\n=== STAGE 2: MERGING AND COMPREHENSIVE QC ===\n")

# merge tmas
merged_obj <- Reduce(function(x,y) merge(x,y), tma_list)

cat("Merged object dimensions:", dim(merged_obj), "\n")
cat("Total spots:", ncol(merged_obj), "\n")
cat("Genotype distribution:\n")
print(table(merged_obj$genotype, useNA = "ifany"))

# Recalculate QC metrics on merged object
merged_obj <- calculate_qc_metrics(merged_obj)

# =============================================================================
# STAGE 3: PRE-FILTERING QC ANALYSIS
# =============================================================================

cat("\n=== STAGE 3: PRE-FILTERING QC ANALYSIS ===\n")

DefaultAssay(merged_obj) <- "Spatial.008um"

# QC metrics by genotype
p1 <- VlnPlot(merged_obj, features = c("nFeature_Spatial.008um", "nCount_Spatial.008um"), 
              group.by = "genotype", ncol = 2, pt.size = 0) + 
  plot_annotation(title = "QC Metrics by Genotype - Pre-filtering")

p2 <- VlnPlot(merged_obj, features = c("percent.mt", "percent.ribo"), 
              group.by = "genotype", ncol = 2, pt.size = 0)

print(p1 / p2)

# QC metrics by TMA (batch effects)
p3 <- VlnPlot(merged_obj, features = c("nFeature_Spatial.008um", "nCount_Spatial.008um"), 
              group.by = "TMA", ncol = 2, pt.size = 0) + 
  plot_annotation(title = "QC Metrics by TMA - Batch Effects Check")

p4 <- VlnPlot(merged_obj, features = c("percent.mt", "percent.ribo"), 
              group.by = "TMA", ncol = 2, pt.size = 0)

print(p3 / p4)

# Endocrine marker validation in Loupe annotations
if(sum(!is.na(merged_obj$Loupe_Cluster)) > 0) {
  cat("Loupe cluster distribution:\n")
  print(table(merged_obj$Loupe_Cluster, useNA = "ifany"))
  
  # Check endocrine markers in annotated regions vs rest
  merged_obj$is_endocrine_annotated <- !is.na(merged_obj$Loupe_Cluster)
  
  endocrine_markers <- c("Ins1", "Ins2", "Gcg", "Sst", "Ppy", "Chga", "Chgb")
  available_markers <- intersect(endocrine_markers, rownames(merged_obj))
  
  if(length(available_markers) > 0) {
    p5 <- VlnPlot(merged_obj, features = available_markers[1:min(4, length(available_markers))], 
                  group.by = "is_endocrine_annotated", ncol = 2, pt.size = 0) +
      plot_annotation(title = "Endocrine Markers: Annotated vs Non-annotated Regions")
    print(p5)
  }
}

# =============================================================================
# STAGE 4: FILTERING DECISIONS AND APPLICATION
# =============================================================================

cat("\n=== STAGE 4: FILTERING ===\n")

# Print current distributions to inform filtering decisions
cat("Current QC metric distributions:\n")
cat("nCount_Spatial.008um range:", range(merged_obj$nCount_Spatial.008um), "\n")
cat("nFeature_Spatial.008um range:", range(merged_obj$nFeature_Spatial.008um), "\n")
cat("percent.mt range:", range(merged_obj$percent.mt), "\n")

# Calculate filtering thresholds
ncount_lower <- quantile(merged_obj$nCount_Spatial.008um, 0.01)
ncount_upper <- quantile(merged_obj$nCount_Spatial.008um, 0.99)
nfeature_lower <- 200
mt_upper <- 20

cat("Filtering thresholds:\n")
cat("nCount range:", ncount_lower, "-", ncount_upper, "\n")
cat("nFeature minimum:", nfeature_lower, "\n")
cat("Mitochondrial maximum:", mt_upper, "%\n")

# Apply spot filtering
spots_before <- ncol(merged_obj)
merged_obj <- subset(merged_obj, 
                    subset = nCount_Spatial.008um > ncount_lower & 
                            nCount_Spatial.008um < ncount_upper &
                            nFeature_Spatial.008um > nfeature_lower & 
                            percent.mt < mt_upper)

spots_after <- ncol(merged_obj)
cat("Spots before filtering:", spots_before, "\n")
cat("Spots after filtering:", spots_after, "\n")
cat("Spots removed:", spots_before - spots_after, "(", round((spots_before - spots_after)/spots_before*100, 1), "%)\n")

# Gene filtering - remove genes expressed in very few spots
DefaultAssay(merged_obj) <- "Spatial.008um"
genes_before <- nrow(merged_obj)

# Keep genes expressed in at least 5 spots
min_spots <- 5
gene_counts <- Matrix::rowSums(GetAssayData(merged_obj, slot = "counts") > 0)
genes_to_keep <- names(gene_counts)[gene_counts >= min_spots]

# Always keep endocrine markers regardless of detection frequency
endocrine_markers <- c("Ins1", "Ins2", "Gcg", "Sst", "Ppy", "Chga", "Chgb", "Isl1", "Neurod1")
endocrine_markers_present <- intersect(endocrine_markers, rownames(merged_obj))
genes_to_keep <- unique(c(genes_to_keep, endocrine_markers_present))

merged_obj <- subset(merged_obj, features = genes_to_keep)

genes_after <- nrow(merged_obj)
cat("Genes before filtering:", genes_before, "\n")
cat("Genes after filtering:", genes_after, "\n")
cat("Genes removed:", genes_before - genes_after, "\n")

# =============================================================================
# STAGE 5: POST-FILTERING VALIDATION
# =============================================================================

cat("\n=== STAGE 5: POST-FILTERING VALIDATION ===\n")

# Check genotype balance after filtering
cat("Post-filtering genotype distribution:\n")
print(table(merged_obj$genotype, useNA = "ifany"))

# Post-filtering QC plots
p6 <- VlnPlot(merged_obj, features = c("nFeature_Spatial.008um", "nCount_Spatial.008um"), 
              group.by = "genotype", ncol = 2, pt.size = 0) + 
  plot_annotation(title = "QC Metrics by Genotype - Post-filtering")

p7 <- VlnPlot(merged_obj, features = c("percent.mt", "percent.ribo"), 
              group.by = "genotype", ncol = 2, pt.size = 0)

print(p6 / p7)

# Spatial QC visualization post-filtering
DefaultAssay(merged_obj) <- "Spatial.008um"
p8 <- SpatialFeaturePlot(merged_obj, features = "nCount_Spatial.008um") + 
  ggtitle("8um nCount Spatial - Post-filtering")
print(p8)

# =============================================================================
# STAGE 6: NORMALIZATION AND ENDOCRINE CELL ANALYSIS
# =============================================================================

cat("\n=== STAGE 6: NORMALIZATION ===\n")

# normalise both assays
DefaultAssay(merged_obj) <- "Spatial.008um"
merged_obj <- NormalizeData(merged_obj)

DefaultAssay(merged_obj) <- "Spatial.016um"
merged_obj <- NormalizeData(merged_obj)

# Focus on endocrine regions if Loupe annotations available
if(sum(!is.na(merged_obj$Loupe_Cluster)) > 0) {
  cat("Subsetting to endocrine-annotated regions\n")
  endocrine_spots <- Cells(merged_obj)[!is.na(merged_obj$Loupe_Cluster)]
  cat("Endocrine spots identified:", length(endocrine_spots), "\n")
  
  # Subset to endocrine regions for main analysis
  merged_obj_endocrine <- subset(merged_obj, cells = endocrine_spots)
  cat("Endocrine subset dimensions:", dim(merged_obj_endocrine), "\n")
  
  # Use endocrine subset for downstream analysis
  analysis_obj <- merged_obj_endocrine
} else {
  cat("No Loupe annotations found, proceeding with full dataset\n")
  analysis_obj <- merged_obj
}

# Endocrine marker visualisations
DefaultAssay(analysis_obj) <- "Spatial.016um"
print(SpatialFeaturePlot(analysis_obj, features = c("Gcg", "Arx"), ncol = 2) + 
      plot_annotation(title = "Alpha cell markers"))
print(SpatialFeaturePlot(analysis_obj, features = c("Ins1", "Ins2"), ncol = 2) + 
      plot_annotation(title = "Beta cell markers"))
print(SpatialFeaturePlot(analysis_obj, features = c("Sst", "Hhex"), ncol = 2) + 
      plot_annotation(title = "Delta cell markers"))
print(SpatialFeaturePlot(analysis_obj, features = c("Ppy", "Chgb"), ncol = 2) + 
      plot_annotation(title = "PP cell markers"))

# =============================================================================
# CLUSTERING AND CELL TYPE IDENTIFICATION
# =============================================================================

cat("\n=== CLUSTERING ANALYSIS ===\n")

# sketched clustering on 8um data
DefaultAssay(analysis_obj) <- "Spatial.008um"
analysis_obj <- FindVariableFeatures(analysis_obj)
analysis_obj <- ScaleData(analysis_obj)
analysis_obj <- SketchData(analysis_obj, ncells = min(50000, ncol(analysis_obj)), 
                          method = "LeverageScore", sketched.assay = "sketch")

DefaultAssay(analysis_obj) <- "sketch"
analysis_obj <- FindVariableFeatures(analysis_obj)
analysis_obj <- ScaleData(analysis_obj)
analysis_obj <- RunPCA(analysis_obj, reduction.name = "pca.sketch")
analysis_obj <- FindNeighbors(analysis_obj, reduction = "pca.sketch", dims = 1:50)
analysis_obj <- FindClusters(analysis_obj, resolution = 3, cluster.name = "seurat_cluster.sketched")
analysis_obj <- RunUMAP(analysis_obj, reduction = "pca.sketch", reduction.name = "umap.sketch", dims = 1:50)

# project clustering onto full data
analysis_obj <- ProjectData(
  object = analysis_obj,
  assay = "Spatial.008um",
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.sketch",
  umap.model = "umap.sketch",
  dims = 1:50,
  refdata = list(seurat_cluster.projected = "seurat_cluster.sketched")
)

# Clustering plots
DefaultAssay(analysis_obj) <- "sketch"
Idents(analysis_obj) <- "seurat_cluster.sketched"
p1 <- DimPlot(analysis_obj, reduction = "umap.sketch", group.by = "seurat_cluster.sketched") + 
  ggtitle("Sketched clustering")

p2 <- DimPlot(analysis_obj, reduction = "umap.sketch", group.by = "genotype") + 
  ggtitle("Genotype distribution")

print(p1 | p2)

DefaultAssay(analysis_obj) <- "Spatial.008um"
Idents(analysis_obj) <- "seurat_cluster.projected"
p3 <- DimPlot(analysis_obj, reduction = "full.umap.sketch") + ggtitle("Projected clustering")
print(p3)

# spatial clusters
print(SpatialDimPlot(analysis_obj, label = TRUE, repel = TRUE, label.size = 4) +
        ggtitle("Spatial clusters - Endocrine regions"))

# =============================================================================
# DIFFERENTIAL EXPRESSION ANALYSIS
# =============================================================================

cat("\n=== DIFFERENTIAL EXPRESSION ANALYSIS ===\n")

# Check cell numbers per genotype
cat("Cells per genotype in analysis:\n")
print(table(analysis_obj$genotype))

# Ensure sufficient cells for comparison
min_cells_per_group <- 50
genotype_counts <- table(analysis_obj$genotype)
if(any(genotype_counts < min_cells_per_group)) {
  warning("Some genotype groups have fewer than", min_cells_per_group, "cells")
}

# DEG analysis comparing genotypes
Idents(analysis_obj) <- analysis_obj$genotype
deg_results <- FindMarkers(analysis_obj, 
                          ident.1 = "Polg_mut", 
                          ident.2 = "WT", 
                          assay = "Spatial.008um", 
                          min.pct = 0.1, 
                          logfc.threshold = 0.25,
                          test.use = "wilcox")

deg_results$gene <- rownames(deg_results)

cat("Total DEGs found:", nrow(deg_results), "\n")
cat("Significant DEGs (padj < 0.05):", sum(deg_results$p_val_adj < 0.05), "\n")

# volcano plot function
plot_volcano <- function(deg_df) {
  deg_df$significant <- deg_df$p_val_adj < 0.05 & abs(deg_df$avg_log2FC) > 0.25
  
  ggplot(deg_df, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
    geom_point(aes(color = significant), alpha=0.6) +
    scale_color_manual(values = c("grey", "red")) +
    theme_minimal() +
    labs(title = "Volcano plot of DEGs (Polg_mut vs WT)", 
         x = "log2 fold change", 
         y = "-log10 adjusted p-value") +
    theme(legend.position = "none") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "red")
}

print(plot_volcano(deg_results))

# =============================================================================
# PATHWAY ENRICHMENT ANALYSIS
# =============================================================================

cat("\n=== PATHWAY ENRICHMENT ANALYSIS ===\n")

# gene conversion for enrichment (gene symbols to entrezid)
gene_list <- deg_results %>%
  filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  pull(gene)

cat("Significant genes for enrichment:", length(gene_list), "\n")

if(length(gene_list) > 0) {
  entrez_ids <- mapIds(org.Mm.eg.db, keys=gene_list, column="ENTREZID", keytype="SYMBOL", multiVals="first")
  entrez_ids <- na.omit(entrez_ids)
  cat("Genes with Entrez IDs:", length(entrez_ids), "\n")
  
  # kegg enrichment
  if(length(entrez_ids) >= 5) {
    kegg_enrich <- enrichKEGG(gene = entrez_ids, organism = 'mmu', pvalueCutoff = 0.05)
    
    if(nrow(kegg_enrich) > 0) {
      print(dotplot(kegg_enrich, showCategory = 15) + ggtitle("KEGG pathway enrichment"))
    } else {
      cat("No significant KEGG pathways found\n")
    }
    
    # gsea analysis
    gene_list_for_gsea <- deg_results$avg_log2FC
    names(gene_list_for_gsea) <- mapIds(org.Mm.eg.db, keys=deg_results$gene, column="ENTREZID", keytype="SYMBOL", multiVals="first")
    gene_list_for_gsea <- gene_list_for_gsea[!is.na(names(gene_list_for_gsea))]
    gene_list_for_gsea <- sort(gene_list_for_gsea, decreasing = TRUE)
    
    msigdbr_sets <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:KEGG")
    gene_sets <- split(msigdbr_sets$gene_symbol, msigdbr_sets$gs_name)
    
    gsea_results <- GSEA(geneList = gene_list_for_gsea,
                         TERM2GENE = data.frame(term=rep(names(gene_sets), lengths(gene_sets)), 
                                               gene=unlist(gene_sets)),
                         pvalueCutoff = 0.05)
    
    if (nrow(gsea_results) > 0) {
      top_pathway <- gsea_results$ID[1]
      print(gseaplot2(gsea_results, geneSetID = top_pathway))
    } else {
      cat("No significant GSEA pathways found.\n")
    }
  } else {
    cat("Too few genes for enrichment analysis\n")
  }
} else {
  cat("No significant DEGs found for enrichment analysis\n")
}

# =============================================================================
# SAVE RESULTS
# =============================================================================

cat("\n=== SAVING RESULTS ===\n")

# Create QC summary
qc_summary <- data.frame(
  metric = c("spots_before_filtering", "spots_after_filtering", "genes_before_filtering", 
             "genes_after_filtering", "endocrine_spots", "total_degs", "significant_degs"),
  value = c(spots_before, spots_after, genes_before, genes_after, 
            ifelse(exists("endocrine_spots"), length(endocrine_spots), ncol(analysis_obj)),
            nrow(deg_results), sum(deg_results$p_val_adj < 0.05))
)

write.csv(qc_summary, "QC_summary.csv", row.names = FALSE)
write.csv(deg_results, "DEG_results_genotype_comparison.csv", row.names = FALSE)
saveRDS(merged_obj, file = "combined_TMA_processed.rds")
saveRDS(analysis_obj, file = "endocrine_analysis_object.rds")

cat("Analysis complete! Files saved:\n")
cat("- QC_summary.csv\n")
cat("- DEG_results_genotype_comparison.csv\n") 
cat("- combined_TMA_processed.rds\n")
cat("- endocrine_analysis_object.rds\n")
