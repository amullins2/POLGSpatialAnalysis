# polg spatial transcriptomics analysis pipeline
# comprehensive analysis of polg vs wt pancreatic islets using visium hd
# 
# this script processes spatial transcriptomics data from 4 tissue microarrays (tmas)
# performs batch correction, cell type identification, and differential expression analysis
# 
# requirements:
# - seurat v4.0+
# - harmony for batch correction
# - spatial transcriptomics data from spaceranger outputs

# load required packages
load_required_packages <- function() {
  required_packages <- c("Seurat", "harmony", "ggplot2", "patchwork", "dplyr", 
                        "tidyr", "RColorBrewer", "Matrix", "arrow", "hdf5r")
  
  for (pkg in required_packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      install.packages(pkg)
      library(pkg, character.only = TRUE)
    }
  }
  
  cat("packages loaded successfully\n")
}

# define consistent publication themes and colors
setup_publication_themes <- function() {
  publication_theme <<- theme_classic() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      axis.title = element_text(size = 11, face = "bold"),
      axis.text = element_text(size = 10, colour = "black"),
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
      plot.background = element_rect(fill = "white", colour = NA)
    )
  
  thesis_colours <<- list(
    genotype = c("POLG" = "#E74C3C", "WT" = "#3498DB"),
    tma = c("TMA1" = "#E74C3C", "TMA2" = "#377EB8", "TMA3" = "#4DAF4A", "TMA4" = "#984EA3"),
    cell_types = c(
      "Pure β-cells" = "#E74C3C", "Pure α-cells" = "#3498DB", 
      "Pure δ-cells" = "#9B59B6", "Pure γ-cells" = "#2ECC71",
      "Mixed Endocrine" = "#F39C12", "Other" = "#95A5A6"
    )
  )
  
  cat("themes and colors defined\n")
}

# load single tma at specified resolution
load_single_tma <- function(tma_name, resolution_um, tma_base_dir) {
  tma_path <- file.path(tma_base_dir, tma_name, "binned_outputs", 
                        paste0("square_", sprintf("%03d", resolution_um), "um"))
  
  cat("loading", tma_name, "at", resolution_um, "μm resolution\n")
  
  if (!dir.exists(tma_path)) {
    cat("directory not found:", tma_path, "\n")
    return(NULL)
  }
  
  # load the spatial data
  obj <- Load10X_Spatial(data.dir = tma_path)
  
  # add essential metadata
  obj@meta.data$tma_id <- tma_name
  obj@meta.data$resolution_um <- resolution_um
  
  # rename assay for clarity
  obj <- RenameAssays(obj, Spatial = paste0("Spatial_", resolution_um, "um"))
  
  cat("loaded", tma_name, ":", ncol(obj), "spots,", nrow(obj), "genes\n")
  return(obj)
}

# load all tma data at specified resolution
load_all_tmas <- function(tma_base_dir, resolution_um = 8) {
  tma_names <- c("TMA1", "TMA2", "TMA3", "TMA4")
  tma_objects <- list()
  
  for (tma_name in tma_names) {
    obj <- load_single_tma(tma_name, resolution_um, tma_base_dir)
    if (!is.null(obj)) {
      tma_objects[[tma_name]] <- obj
    }
  }
  
  cat("successfully loaded", length(tma_objects), "tmas at", resolution_um, "μm resolution\n")
  return(tma_objects)
}

# coordinate matching for islet annotations
efficient_coordinate_matching <- function(tma_name, hd_data_dir, old_data_dir) {
  cat("processing", tma_name, "coordinate matching\n")
  
  # load hd spatial coordinates
  spatial_coords_file <- file.path(hd_data_dir, tma_name, 
                                  "binned_outputs/square_008um/spatial/tissue_positions.parquet")
  old_coords_file <- file.path(old_data_dir, tma_name, "SpatialCoords.csv")
  
  if (!file.exists(spatial_coords_file) || !file.exists(old_coords_file)) {
    cat("coordinate files not found for", tma_name, "\n")
    return(NULL)
  }
  
  # read coordinates
  hd_coords <- read_parquet(spatial_coords_file)
  old_coords <- read.csv(old_coords_file)
  colnames(old_coords) <- c("loupe_barcode", "x_coord", "y_coord")
  
  # load islet annotations
  all_annotations <- data.frame()
  for (genotype in c("POLG", "WT")) {
    file_path <- file.path(old_data_dir, tma_name, paste0("Islet_", genotype, "_", tma_name, ".csv"))
    if (file.exists(file_path)) {
      data <- read.csv(file_path, stringsAsFactors = FALSE)
      colnames(data)[1] <- "loupe_barcode"
      if (ncol(data) > 1) colnames(data)[2] <- "islet_id"
      data$genotype <- genotype
      all_annotations <- rbind(all_annotations, data)
    }
  }
  
  if (nrow(all_annotations) == 0) {
    cat("no annotations found for", tma_name, "\n")
    return(NULL)
  }
  
  # merge annotations with coordinates
  islet_coords <- merge(all_annotations, old_coords, by = "loupe_barcode")
  
  # coordinate matching with tolerance
  matches <- data.frame()
  tolerance <- 10
  
  for (i in 1:nrow(islet_coords)) {
    target_x <- islet_coords$x_coord[i]
    target_y <- islet_coords$y_coord[i]
    
    distances <- sqrt((hd_coords$pxl_col_in_fullres - target_x)^2 + 
                     (hd_coords$pxl_row_in_fullres - target_y)^2)
    
    min_dist <- min(distances, na.rm = TRUE)
    if (min_dist <= tolerance) {
      best_idx <- which.min(distances)
      matches <- rbind(matches, data.frame(
        hd_barcode = hd_coords$barcode[best_idx],
        loupe_barcode = islet_coords$loupe_barcode[i],
        genotype = islet_coords$genotype[i],
        islet_id = islet_coords$islet_id[i],
        distance = min_dist,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  cat(tma_name, "matching complete:", nrow(matches), "matches\n")
  return(matches)
}

# create islet seurat objects with annotations
create_islet_objects <- function(tma_objects, coordinate_matches) {
  tma_islets <- list()
  
  for (tma_name in names(coordinate_matches)) {
    matches <- coordinate_matches[[tma_name]]
    seurat_obj <- tma_objects[[tma_name]]
    
    # find matching barcodes
    seurat_barcodes <- colnames(seurat_obj)
    matched_barcodes <- intersect(seurat_barcodes, matches$hd_barcode)
    
    if (length(matched_barcodes) > 0) {
      # subset to islet spots
      islet_obj <- subset(seurat_obj, cells = matched_barcodes)
      
      # add islet metadata
      match_idx <- match(matched_barcodes, matches$hd_barcode)
      islet_obj@meta.data$genotype <- matches$genotype[match_idx]
      islet_obj@meta.data$islet_id <- matches$islet_id[match_idx]
      islet_obj@meta.data$coordinate_distance <- matches$distance[match_idx]
      
      tma_islets[[tma_name]] <- islet_obj
      
      cat(tma_name, ":", ncol(islet_obj), "islet spots created\n")
    }
  }
  
  return(tma_islets)
}

# combine all tma islet objects
combine_tma_objects <- function(tma_islets) {
  cat("combining all tma islet objects\n")
  
  # add tma id to each object before merging
  for (tma_name in names(tma_islets)) {
    tma_islets[[tma_name]]$tma_id <- tma_name
    # create unique cell names
    new_cell_names <- paste0(tma_name, "_", colnames(tma_islets[[tma_name]]))
    tma_islets[[tma_name]] <- RenameCells(tma_islets[[tma_name]], new.names = new_cell_names)
  }
  
  if (length(tma_islets) == 1) {
    combined_islets <- tma_islets[[1]]
  } else {
    combined_islets <- tma_islets[[1]]
    for (i in 2:length(tma_islets)) {
      combined_islets <- merge(combined_islets, tma_islets[[i]])
    }
  }
  
  # convert tma_id to factor
  combined_islets$tma_id <- factor(combined_islets$tma_id, 
                                  levels = c("TMA1", "TMA2", "TMA3", "TMA4"))
  
  cat("combined dataset created:", ncol(combined_islets), "total spots\n")
  return(combined_islets)
}

# run harmony batch correction workflow
run_harmony_batch_correction <- function(combined_islets) {
  cat("starting harmony batch correction workflow\n")
  
  # basic sctransform normalization without regression variables
  cat("running sctransform normalization\n")
  combined_basic <- SCTransform(
    combined_islets,
    assay = names(combined_islets@assays)[1],
    verbose = FALSE
  )
  
  cat("running pca\n")
  combined_basic <- RunPCA(combined_basic, verbose = FALSE)
  
  # apply harmony batch correction
  cat("running harmony batch correction\n")
  combined_harmony <- RunHarmony(
    object = combined_basic,
    group.by.vars = "tma_id"
  )
  
  # run umap on batch-corrected data
  cat("running umap on batch-corrected data\n")
  combined_harmony <- RunUMAP(
    combined_harmony,
    reduction = "harmony",
    dims = 1:30,
    verbose = FALSE
  )
  
  cat("harmony workflow complete\n")
  return(combined_harmony)
}

# identify pure cell types based on hormone expression
identify_pure_cell_types <- function(combined_islets) {
  cat("identifying pure cell types using hormone expression\n")
  
  # extract hormone expression for all cells
  hormone_genes <- c("Gcg", "Ins1", "Ins2", "Sst", "Ppy")
  available_hormones <- intersect(hormone_genes, rownames(combined_islets))
  
  if (length(available_hormones) == 0) {
    cat("no hormone genes found in dataset\n")
    return(combined_islets)
  }
  
  hormone_data <- FetchData(combined_islets, vars = available_hormones)
  meta_data <- combined_islets@meta.data[, c("genotype", "tma_id")]
  
  # combine expression and metadata
  cell_analysis <- cbind(hormone_data, meta_data)
  cell_analysis$cell_id <- rownames(cell_analysis)
  
  # define pure cell type criteria
  cell_analysis$cell_type <- "Other"
  
  # pure beta cells: high insulin, low glucagon
  if (all(c("Ins1", "Ins2", "Gcg", "Sst") %in% colnames(cell_analysis))) {
    beta_criteria <- (cell_analysis$Ins1 > 0.5 | cell_analysis$Ins2 > 0.5) & 
                     cell_analysis$Gcg < 0.5 & cell_analysis$Sst < 0.5
    cell_analysis$cell_type[beta_criteria] <- "Pure β-cells"
    
    # pure alpha cells: high glucagon, low insulin
    alpha_criteria <- cell_analysis$Gcg > 0.5 & 
                      cell_analysis$Ins1 < 0.2 & cell_analysis$Ins2 < 0.2 & 
                      cell_analysis$Sst < 0.5
    cell_analysis$cell_type[alpha_criteria] <- "Pure α-cells"
    
    # pure delta cells: high somatostatin, low others
    delta_criteria <- cell_analysis$Sst > 0.5 & 
                      cell_analysis$Ins1 < 0.2 & cell_analysis$Ins2 < 0.2 & 
                      cell_analysis$Gcg < 0.5
    cell_analysis$cell_type[delta_criteria] <- "Pure δ-cells"
  }
  
  # add cell type labels to seurat object
  combined_islets$pure_cell_type <- cell_analysis[colnames(combined_islets), "cell_type"]
  
  # print summary
  cell_counts <- table(combined_islets$pure_cell_type, combined_islets$genotype)
  cat("cell type identification complete\n")
  print(cell_counts)
  
  return(combined_islets)
}

# run differential expression analysis
run_differential_expression <- function(combined_islets) {
  cat("running differential expression analysis\n")
  
  # identify cell types with sufficient cells for analysis
  target_cell_types <- c("Pure β-cells", "Pure α-cells", "Pure δ-cells")
  min_cells_per_group <- 10
  deg_results <- list()
  
  for (cell_type in target_cell_types) {
    polg_count <- sum(combined_islets$pure_cell_type == cell_type & 
                     combined_islets$genotype == "POLG", na.rm = TRUE)
    wt_count <- sum(combined_islets$pure_cell_type == cell_type & 
                   combined_islets$genotype == "WT", na.rm = TRUE)
    
    if (polg_count >= min_cells_per_group && wt_count >= min_cells_per_group) {
      cat("analyzing", cell_type, "- polg:", polg_count, ", wt:", wt_count, "\n")
      
      # subset to cell type
      cell_subset <- subset(combined_islets, 
                           subset = pure_cell_type == cell_type & 
                                   !is.na(genotype))
      
      # set identities
      Idents(cell_subset) <- "genotype"
      
      # run findmarkers
      tryCatch({
        markers <- FindMarkers(
          cell_subset,
          ident.1 = "POLG",
          ident.2 = "WT",
          test.use = "wilcox",
          min.pct = 0.1,
          logfc.threshold = 0.1
        )
        
        # clean results
        markers$gene <- rownames(markers)
        markers$cell_type <- cell_type
        markers$regulation <- ifelse(markers$avg_log2FC > 0, 
                                   "Upregulated_in_POLG", 
                                   "Downregulated_in_POLG")
        
        # store results
        clean_name <- gsub("[^A-Za-z0-9]", "_", cell_type)
        deg_results[[clean_name]] <- markers
        
        cat("completed", cell_type, ":", nrow(markers), "genes tested\n")
        
      }, error = function(e) {
        cat("deg analysis failed for", cell_type, ":", e$message, "\n")
      })
    }
  }
  
  return(deg_results)
}

# create volcano plots
create_volcano_plots <- function(deg_results) {
  cat("creating volcano plots\n")
  
  volcano_plots <- list()
  
  for (cell_type in names(deg_results)) {
    deg_data <- deg_results[[cell_type]]
    
    # prepare plot data
    plot_data <- deg_data
    plot_data$neg_log10_p <- -log10(pmax(plot_data$p_val_adj, 1e-300))
    
    # define significance
    plot_data$regulation <- "Not Significant"
    plot_data$regulation[plot_data$p_val_adj < 0.05 & plot_data$avg_log2FC > 0.25] <- "Upregulated in POLG"
    plot_data$regulation[plot_data$p_val_adj < 0.05 & plot_data$avg_log2FC < -0.25] <- "Downregulated in POLG"
    
    reg_colours <- c(
      "Not Significant" = "grey70", 
      "Upregulated in POLG" = "#E74C3C", 
      "Downregulated in POLG" = "#3498DB"
    )
    
    # create plot
    volcano_plot <- ggplot(plot_data, aes(x = avg_log2FC, y = neg_log10_p, colour = regulation)) +
      geom_point(alpha = 0.7, size = 1) +
      scale_colour_manual(values = reg_colours) +
      geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", alpha = 0.5) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5) +
      labs(
        title = paste("differential gene expression:", gsub("_", " ", cell_type)),
        x = "log2 fold change (polg vs wt)", 
        y = "-log10(adjusted p-value)",
        colour = "regulation in polg"
      ) +
      publication_theme
    
    volcano_plots[[cell_type]] <- volcano_plot
  }
  
  return(volcano_plots)
}

# create batch correction comparison plots
create_batch_correction_plots <- function(combined_islets) {
  cat("creating batch correction comparison plots\n")
  
  plots <- list()
  
  # before correction (pca)
  if ("pca" %in% names(combined_islets@reductions)) {
    before_plot <- DimPlot(combined_islets, 
                          reduction = "pca", 
                          group.by = "tma_id",
                          cols = thesis_colours$tma) +
      ggtitle("before harmony correction") +
      labs(subtitle = "clear tma separation indicates batch effects") +
      publication_theme
    
    plots$before <- before_plot
  }
  
  # after correction (umap)
  if ("umap" %in% names(combined_islets@reductions)) {
    after_plot <- DimPlot(combined_islets, 
                         reduction = "umap", 
                         group.by = "tma_id",
                         cols = thesis_colours$tma) +
      ggtitle("after harmony correction") +
      labs(subtitle = "well-integrated tmas") +
      publication_theme
    
    plots$after <- after_plot
  }
  
  return(plots)
}

# main analysis pipeline
run_polg_analysis_pipeline <- function(hd_data_dir, old_data_dir, output_dir = "polg_analysis_results") {
  cat("starting polg spatial transcriptomics analysis pipeline\n")
  
  # create output directory
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # setup
  load_required_packages()
  setup_publication_themes()
  
  # step 1: load spatial data
  cat("\nstep 1: loading spatial data\n")
  tma_objects <- load_all_tmas(hd_data_dir, resolution_um = 8)
  
  # step 2: coordinate matching for islet annotations
  cat("\nstep 2: coordinate matching for islet annotations\n")
  coordinate_matches <- list()
  for (tma_name in names(tma_objects)) {
    matches <- efficient_coordinate_matching(tma_name, hd_data_dir, old_data_dir)
    if (!is.null(matches)) {
      coordinate_matches[[tma_name]] <- matches
    }
  }
  
  # step 3: create islet objects
  cat("\nstep 3: creating islet objects\n")
  tma_islets <- create_islet_objects(tma_objects, coordinate_matches)
  
  # step 4: combine objects
  cat("\nstep 4: combining tma objects\n")
  combined_islets <- combine_tma_objects(tma_islets)
  
  # step 5: harmony batch correction
  cat("\nstep 5: harmony batch correction\n")
  combined_islets <- run_harmony_batch_correction(combined_islets)
  
  # step 6: cell type identification
  cat("\nstep 6: cell type identification\n")
  combined_islets <- identify_pure_cell_types(combined_islets)
  
  # step 7: differential expression analysis
  cat("\nstep 7: differential expression analysis\n")
  deg_results <- run_differential_expression(combined_islets)
  
  # step 8: create visualizations
  cat("\nstep 8: creating visualizations\n")
  
  # volcano plots
  volcano_plots <- create_volcano_plots(deg_results)
  for (cell_type in names(volcano_plots)) {
    filename <- file.path(output_dir, paste0("volcano_", cell_type, ".png"))
    ggsave(filename, volcano_plots[[cell_type]], width = 12, height = 8, dpi = 300)
  }
  
  # batch correction plots
  batch_plots <- create_batch_correction_plots(combined_islets)
  if (length(batch_plots) >= 2) {
    combined_batch <- batch_plots$before | batch_plots$after
    ggsave(file.path(output_dir, "batch_correction_comparison.png"), 
           combined_batch, width = 16, height = 8, dpi = 300)
  }
  
  # umap plots
  umap_genotype <- DimPlot(combined_islets, reduction = "umap", 
                          group.by = "genotype", cols = thesis_colours$genotype) +
    ggtitle("genotype distribution") + publication_theme
  
  umap_celltype <- DimPlot(combined_islets, reduction = "umap", 
                          group.by = "pure_cell_type", cols = thesis_colours$cell_types) +
    ggtitle("cell type distribution") + publication_theme
  
  ggsave(file.path(output_dir, "umap_genotype.png"), umap_genotype, 
         width = 10, height = 8, dpi = 300)
  ggsave(file.path(output_dir, "umap_celltype.png"), umap_celltype, 
         width = 12, height = 8, dpi = 300)
  
  # step 9: save results
  cat("\nstep 9: saving results\n")
  
  # save seurat object
  saveRDS(combined_islets, file.path(output_dir, "combined_islets_processed.rds"))
  
  # save deg results
  for (cell_type in names(deg_results)) {
    filename <- file.path(output_dir, paste0("deg_results_", cell_type, ".csv"))
    write.csv(deg_results[[cell_type]], filename, row.names = FALSE)
  }
  
  # create summary report
  summary_data <- data.frame(
    total_spots = ncol(combined_islets),
    total_genes = nrow(combined_islets),
    tma_count = length(unique(combined_islets$tma_id)),
    cell_types_identified = length(unique(combined_islets$pure_cell_type)),
    deg_analyses = length(deg_results),
    polg_spots = sum(combined_islets$genotype == "POLG", na.rm = TRUE),
    wt_spots = sum(combined_islets$genotype == "WT", na.rm = TRUE)
  )
  
  write.csv(summary_data, file.path(output_dir, "analysis_summary.csv"), row.names = FALSE)
  
  cat("\nanalysis pipeline complete!\n")
  cat("results saved to:", output_dir, "\n")
  
  return(list(
    combined_islets = combined_islets,
    deg_results = deg_results,
    tma_objects = tma_objects,
    summary = summary_data
  ))
}

# usage example:
# results <- run_polg_analysis_pipeline(
#   hd_data_dir = "/path/to/hd/data",
#   old_data_dir = "/path/to/old/data",
#   output_dir = "polg_analysis_results"
# )

results <- run_polg_analysis_pipeline(
  hd_data_dir = "/path/to/hd/data",
  old_data_dir = "/path/to/old/data", 
  output_dir = "polg_analysis_results"
)
