# spatial transcriptomics analysis pipeline for mouse pancreatic tmas
#
# this script performs the following tasks:
# - installs and loads required cran and bioconductor packages
# - loads spatial transcriptomics data from multiple tmas with loupe cluster annotations
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
# adjust file paths for metadata, tma data, loupe csvs, and gene lists accordingly.
#
# author: alana mullins
# date: 18/05/2025
# =============================================================================

# install cran packages, if not there
cran_packages <- c("Seurat", "dplyr", "ggplot2", "patchwork")
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
library(clusterProfiler)
library(org.Mm.eg.db)   
library(msigdbr)        
library(fgsea)          

# load metadata
metadata <- read.csv("/sample_metadata.csv", stringsAsFactors = FALSE)

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
  object <- Load10X_Spatial(data.dir = tma_dir, bin.size = c(8, 16))
  
  # loupe clusters
  loupe_annot <- read.csv(loupe_csv_path)
  rownames(loupe_annot) <- loupe_annot$barcode
  common_barcodes <- intersect(Cells(object), rownames(loupe_annot))
  object$Loupe_Cluster <- NA
  object$Loupe_Cluster[common_barcodes] <- loupe_annot[common_barcodes, "cluster"]
  
  # extract core_id & add metadata
  object$core_id <- sapply(strsplit(Cells(object), "-"), `[`, 1)
  for (field in c("mouse_id", "genotype", "sex", "condition")) {
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

# merge tmas
merged_obj <- Reduce(function(x,y) merge(x,y), tma_list)

# qc before normalisation
DefaultAssay(merged_obj) <- "Spatial.008um"
p1 <- VlnPlot(merged_obj, features = "nCount_Spatial.008um", pt.size = 0) + NoLegend() + ggtitle("8um nCount")
p2 <- SpatialFeaturePlot(merged_obj, features = "nCount_Spatial.008um") + ggtitle("8um nCount Spatial")
print(p1 | p2)

# normalise
DefaultAssay(merged_obj) <- "Spatial.008um"
merged_obj <- NormalizeData(merged_obj)

DefaultAssay(merged_obj) <- "Spatial.016um"
merged_obj <- NormalizeData(merged_obj)

# marker visualisations
DefaultAssay(merged_obj) <- "Spatial.016um"
print(SpatialFeaturePlot(merged_obj, features = c("Gcg", "Arx", "Ttr")) + ggtitle("Alpha markers"))
print(SpatialFeaturePlot(merged_obj, features = c("Ins1", "Pdx1", "Nkx6-1")) + ggtitle("Beta markers"))
print(SpatialFeaturePlot(merged_obj, features = c("Sst", "Hhex", "Rbp4")) + ggtitle("Delta markers"))
print(SpatialFeaturePlot(merged_obj, features = c("Ppy", "Pax6", "Chgb")) + ggtitle("PP markers"))

# sketched clustering
DefaultAssay(merged_obj) <- "Spatial.008um"
merged_obj <- FindVariableFeatures(merged_obj)
merged_obj <- ScaleData(merged_obj)
merged_obj <- SketchData(merged_obj, ncells = 50000, method = "LeverageScore", sketched.assay = "sketch")

DefaultAssay(merged_obj) <- "sketch"
merged_obj <- FindVariableFeatures(merged_obj)
merged_obj <- ScaleData(merged_obj)
merged_obj <- RunPCA(merged_obj, reduction.name = "pca.sketch")
merged_obj <- FindNeighbors(merged_obj, reduction = "pca.sketch", dims = 1:50)
merged_obj <- FindClusters(merged_obj, resolution = 3, cluster.name = "seurat_cluster.sketched")
merged_obj <- RunUMAP(merged_obj, reduction = "pca.sketch", reduction.name = "umap.sketch", dims = 1:50)

# project clustering onto full data
merged_obj <- ProjectData(
  object = merged_obj,
  assay = "Spatial.008um",
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.sketch",
  umap.model = "umap.sketch",
  dims = 1:50,
  refdata = list(seurat_cluster.projected = "seurat_cluster.sketched")
)

# umap plots
DefaultAssay(merged_obj) <- "sketch"
Idents(merged_obj) <- "seurat_cluster.sketched"
p1 <- DimPlot(merged_obj, reduction = "umap.sketch") + ggtitle("Sketched clustering")

DefaultAssay(merged_obj) <- "Spatial.008um"
Idents(merged_obj) <- "seurat_cluster.projected"
p2 <- DimPlot(merged_obj, reduction = "full.umap.sketch") + ggtitle("Projected clustering")

print(p1 | p2)

# spatial clusters
print(SpatialDimPlot(merged_obj, label = TRUE, repel = TRUE, label.size = 4) +
        ggtitle("Merged tmas - spatial clusters"))

# highlight endocrine clusters (change cluster ids accordingly)
endocrine_clusters <- c(0,4,32,34,35)
Idents(merged_obj) <- "seurat_cluster.projected"
cells_endocrine <- WhichCells(merged_obj, idents = endocrine_clusters)

p_highlight <- SpatialDimPlot(merged_obj,
                              cells.highlight = unlist(cells_endocrine),
                              cols.highlight = c("#FFFF00", "grey50"),
                              facet.highlight = TRUE) + NoLegend()
print(p_highlight)

# subset endocrine cells for deg analysis
endocrine_obj <- subset(merged_obj, idents = endocrine_clusters)

# deg analysis comparing genotype 

Idents(endocrine_obj) <- endocrine_obj$genotype
deg_results <- FindMarkers(endocrine_obj, ident.1 = "Polg_mut", ident.2 = "WT", assay = "Spatial.008um", min.pct = 0.1, logfc.threshold = 0.25)

deg_results$gene <- rownames(deg_results)

# volcano plot function
plot_volcano <- function(deg_df) {
  ggplot(deg_df, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
    geom_point(aes(color = p_val_adj < 0.05 & abs(avg_log2FC) > 0.25), alpha=0.6) +
    scale_color_manual(values = c("grey", "red")) +
    theme_minimal() +
    labs(title = "Volcano plot of degs", x = "log2 fold change", y = "-log10 adjusted p-value") +
    theme(legend.position = "none")
}

print(plot_volcano(deg_results))

# gene conversion for enrichment (gene symbols to entrezid)

gene_list <- deg_results %>%
  filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  pull(gene)

entrez_ids <- mapIds(org.Mm.eg.db, keys=gene_list, column="ENTREZID", keytype="SYMBOL", multiVals="first")

entrez_ids <- na.omit(entrez_ids)

# mitocarta enrichment 

# load mitocarta gene list
mitocarta_genes <- read.csv("/mitocarta_genes.csv", header = TRUE, stringsAsFactors = FALSE)$gene_symbol

mito_gene_entrez <- mapIds(org.Mm.eg.db, keys=mitocarta_genes, column="ENTREZID", keytype="SYMBOL", multiVals="first")
mito_gene_entrez <- na.omit(mito_gene_entrez)

# perform enrichment (hypergeometric test)
mitocarta_enrich <- enricher(entrez_ids, TERM2GENE = data.frame(term = "MitoCarta", gene = mito_gene_entrez))

# bar plot for mitocarta enrichment
plot_mitocarta_enrichment <- function(enrich_result) {
  df <- as.data.frame(enrich_result)
  df <- df[order(df$p.adjust), ][1:min(10, nrow(df)), ]
  
  ggplot(df, aes(x = reorder(Description, -p.adjust), y = -log10(p.adjust))) +
    geom_bar(stat="identity", fill="steelblue") +
    coord_flip() +
    labs(title = "Top MitoCarta enriched terms", x = "", y = "-log10 adjusted p-value") +
    theme_minimal()
}

print(plot_mitocarta_enrichment(mitocarta_enrich))

# kegg enrichment

kegg_enrich <- enrichKEGG(gene = entrez_ids, organism = 'mmu', pvalueCutoff = 0.05)

print(dotplot(kegg_enrich, showCategory = 15) + ggtitle("KEGG pathway enrichment"))
print(emapplot(kegg_enrich))

# gsea analysis

# prepare ranked gene list for gsea (named vector, decreasing order)
gene_list_for_gsea <- deg_results$avg_log2FC
names(gene_list_for_gsea) <- mapIds(org.Mm.eg.db, keys=deg_results$gene, column="ENTREZID", keytype="SYMBOL", multiVals="first")
gene_list_for_gsea <- gene_list_for_gsea[!is.na(names(gene_list_for_gsea))]
gene_list_for_gsea <- sort(gene_list_for_gsea, decreasing = TRUE)

msigdbr_sets <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:KEGG")

gene_sets <- split(msigdbr_sets$gene_symbol, msigdbr_sets$gs_name)

gsea_results <- GSEA(geneList = gene_list_for_gsea,
                     TERM2GENE = data.frame(term=rep(names(gene_sets), lengths(gene_sets)), gene=unlist(gene_sets)),
                     pvalueCutoff = 0.05)

# plot gsea for top pathway
if (nrow(gsea_results) > 0) {
  top_pathway <- gsea_results$ID[1]
  print(gseaplot2(gsea_results, geneSetID = top_pathway))
} else {
  message("No significant gsea pathways found.")
}

# save combined object and deg results 

saveRDS(merged_obj, file = "combined_TMA_processed.rds")
write.csv(deg_results, "DEG_results_genotype_comparison.csv", row.names = FALSE)
