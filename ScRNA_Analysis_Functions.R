# =============================================================================
# ScRNA-Seq Single-Cell Analysis — Custom Functions Library
# =============================================================================
# Author : Ellie Malcolm et al.
# Date   : 2025-03
# Description: Reusable functions for QC, preprocessing, clustering,
#              pseudobulk generation, annotation, and data export.
# =============================================================================
#
# TABLE OF CONTENTS
# -----------------
#  1. QC and Visualisation
#       load_cellbender_filtered_h5()
#       plot_qc_violin_grid()
#       resumen_nFeature_plot()
#       safe_vln()
#       mostrar_tabla()
#       plot_integrated_clusters()   [sourced from custom_seurat.R]
#
#  2. Preprocessing and Doublet Detection
#       preprocesar_y_doubletfinder()
#       doubletfinder_pipeline()
#       process_sample()
#
#  3. Pseudobulk
#       generate_pseudobulk()
#       normalizar_bulk_pseudobulk()
#       clasificar_residuos()
#
#  4. Clustering Utilities
#       unificar_nombres()
#       unir_layers_counts()
#
#  5. Marker-Based and Reference-Based Annotation
#       find_markers()
#       annotate_by_markers()
#       annotate_by_reference()
#
#  6. Export
#       exportar_para_scanpy()
#
# =============================================================================


# =============================================================================
# 1. QC AND VISUALISATION
# =============================================================================

#' Load CellBender Filtered HDF5 Data
#'
#' Reads a filtered expression matrix from CellBender HDF5 output and returns
#' a Seurat object with raw counts.
#'
#' @param h5_path  Path to the filtered HDF5 file produced by CellBender.
#' @param project  Project name stored in Seurat metadata.
#' @return A Seurat object.
#' @export
load_cellbender_filtered_h5 <- function(h5_path, project = "Sample") {

  f <- H5File$new(h5_path, mode = "r")

  message("Reading CSR components...")
  data    <- f[["matrix/data"]]$read()
  indices <- f[["matrix/indices"]]$read()
  indptr  <- f[["matrix/indptr"]]$read()
  shape   <- f[["matrix/shape"]]$read()

  message("Reading gene IDs and barcodes...")
  gene_ids <- f[["matrix/features/id"]]$read()
  barcodes <- f[["matrix/barcodes"]]$read()

  message("Creating sparse gene x cell matrix...")
  mat <- new("dgCMatrix",
             x        = as.numeric(data),
             i        = indices,
             p        = indptr,
             Dim      = shape,
             Dimnames = list(gene_ids, barcodes))

  seu <- CreateSeuratObject(counts = mat, project = project)
  return(seu)
}


#' QC Violin Plot Grid
#'
#' Generates a 4-panel violin plot (nFeature_RNA, nCount_RNA, percent.mt,
#' percent.cp) for a single sample, grouped by condition label.
#'
#' @param obj1  Seurat object.
#' @param label Condition label shown in the plot title.
#' @param color Fill colour for the violins.
#' @return A ggplot object.
#' @export
plot_qc_violin_grid <- function(obj1, label, color) {

  n1       <- ncol(obj1)
  obj1$cond <- label

  p1 <- VlnPlot(
    obj1,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.cp"),
    pt.size  = 0.1,
    ncol     = 4,
    group.by = "cond",
    cols     = color
  ) +
    ggtitle(paste0(label, " (", n1, " cells)")) +
    theme_minimal(base_size = 12)

  return(p1)
}


#' Summary of nFeature_RNA Distribution
#'
#' Produces a boxplot with overlaid jitter points together with quartile and
#' quintile summary tables for a list of Seurat objects.
#'
#' @param obj_list  Named list of Seurat objects.
#' @param etiquetas Character vector of labels (one per object).
#' @param colores   Named colour vector; defaults to ColorBrewer Set2 palette.
#' @export
resumen_nFeature_plot <- function(obj_list, etiquetas = NULL, colores = NULL) {

  if (is.null(etiquetas)) etiquetas <- paste0("Group", seq_along(obj_list))
  if (length(etiquetas) != length(obj_list)) stop("Labels must match objects.")

  if (is.null(colores)) {
    colores <- c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854")[
      seq_along(obj_list)
    ]
    names(colores) <- etiquetas
  }

  lista_df <- lapply(seq_along(obj_list), function(i) {
    obj <- obj_list[[i]]
    data.frame(nFeature_RNA = obj@meta.data$nFeature_RNA, grupo = etiquetas[i])
  })

  meta_comb        <- bind_rows(lista_df)
  meta_comb$grupo  <- factor(meta_comb$grupo, levels = etiquetas)

  p_box <- ggplot(meta_comb, aes(x = grupo, y = nFeature_RNA, fill = grupo)) +
    geom_boxplot(outlier.shape = NA, width = 0.6) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
    scale_fill_manual(values = colores) +
    labs(title = "nFeature_RNA Distribution", x = "Condition", y = "nFeature_RNA") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none")

  cuartiles <- meta_comb %>%
    group_by(grupo) %>%
    summarise(
      Min    = quantile(nFeature_RNA, 0.00),
      Q1     = quantile(nFeature_RNA, 0.25),
      Median = quantile(nFeature_RNA, 0.50),
      Q3     = quantile(nFeature_RNA, 0.75),
      Max    = quantile(nFeature_RNA, 1.00),
      .groups = "drop"
    ) %>%
    arrange(factor(grupo, levels = etiquetas))

  quintiles <- meta_comb %>%
    group_by(grupo) %>%
    summarise(
      `0%`   = quantile(nFeature_RNA, 0.0),
      `20%`  = quantile(nFeature_RNA, 0.2),
      `40%`  = quantile(nFeature_RNA, 0.4),
      `60%`  = quantile(nFeature_RNA, 0.6),
      `80%`  = quantile(nFeature_RNA, 0.8),
      `100%` = quantile(nFeature_RNA, 1.0),
      .groups = "drop"
    ) %>%
    arrange(factor(grupo, levels = etiquetas))

  tabla_cuartiles <- tableGrob(cuartiles)
  tabla_quintiles <- tableGrob(quintiles)

  panel_tablas <- plot_grid(
    ggdraw() + draw_label("Quartiles", fontface = "bold", size = 13),
    ggdraw() + draw_grob(tabla_cuartiles),
    ggdraw() + draw_label("Quintiles", fontface = "bold", size = 13),
    ggdraw() + draw_grob(tabla_quintiles),
    ncol        = 1,
    rel_heights = c(0.15, 1, 0.15, 1)
  )

  final_plot <- plot_grid(p_box, panel_tablas, ncol = 2, rel_widths = c(1.5, 1))
  print(final_plot)
}


#' Safe Violin Plot for RMarkdown
#'
#' Wrapper around VlnPlot with manual colour scale, designed for clean
#' rendering inside RMarkdown documents.
#'
#' @param obj     Seurat object.
#' @param feature Gene or metadata feature to plot.
#' @param colors  Named colour vector for orig.ident groups.
#' @return A ggplot object.
#' @export
safe_vln <- function(obj, feature, colors) {
  p <- VlnPlot(obj, features = feature, group.by = "orig.ident", pt.size = 0)
  p <- p + scale_fill_manual(values = colors)
  return(p)
}


#' Display Annotation Table as Grid
#'
#' Creates a side-by-side count table comparing cell-type annotations between
#' two vectors (e.g., filtered vs CellBender-corrected) and renders it using
#' grid graphics.
#'
#' @param filtered_vec   Annotation vector for the filtered dataset.
#' @param cellbender_vec Annotation vector for the CellBender dataset.
#' @param titulo         Table title (currently informational only).
#' @export
mostrar_tabla <- function(filtered_vec, cellbender_vec, titulo = "Annotations") {

  t1       <- table(filtered_vec)
  t2       <- table(cellbender_vec)
  all_types <- union(names(t1), names(t2))

  df <- data.frame(
    celltype   = all_types,
    filtered   = as.integer(t1[all_types]),
    cellbender = as.integer(t2[all_types]),
    stringsAsFactors = FALSE
  )
  df[is.na(df)] <- 0

  total_row <- data.frame(
    celltype   = "Total",
    filtered   = sum(df$filtered),
    cellbender = sum(df$cellbender),
    stringsAsFactors = FALSE
  )
  df <- rbind(df, total_row)

  grid.newpage()
  grid.draw(tableGrob(df, rows = NULL, theme = ttheme_minimal()))
}


# =============================================================================
# 2. PREPROCESSING AND DOUBLET DETECTION
# =============================================================================

#' Normalize, Scale, PCA, and Run DoubletFinder
#'
#' Simplified preprocessing pipeline followed by DoubletFinder doublet
#' detection. Returns the Seurat object with doublet classifications in
#' metadata but does NOT filter singlets.
#'
#' @param seurat_obj             Seurat object.
#' @param pcs                    Integer vector of PCs to use (default 1:20).
#' @param expected_doublet_rate  Expected fraction of doublets (default 0.075).
#' @param project_id             Project label (informational).
#' @return Seurat object with DoubletFinder classifications.
#' @export
preprocesar_y_doubletfinder <- function(seurat_obj,
                                        pcs                    = 1:20,
                                        expected_doublet_rate  = 0.075,
                                        project_id             = "sample") {
  message("Normalizing...")
  seurat_obj <- NormalizeData(seurat_obj)

  message("Finding variable features...")
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst",
                                     nfeatures = 2000)
  message("Scaling...")
  seurat_obj <- ScaleData(seurat_obj)

  message("Running PCA...")
  seurat_obj <- RunPCA(seurat_obj, npcs = max(pcs))

  message("Running DoubletFinder...")
  sweep.res   <- paramSweep(seurat_obj, PCs = pcs, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
  bcmvn       <- find.pK(sweep.stats)
  best.pK     <- as.numeric(as.character(bcmvn[which.max(bcmvn$BCmetric), "pK"]))
  nExp        <- round(expected_doublet_rate * ncol(seurat_obj))

  seurat_obj <- doubletFinder(
    seurat_obj,
    PCs        = pcs,
    pN         = 0.25,
    pK         = best.pK,
    nExp       = nExp,
    reuse.pANN = FALSE,
    sct        = FALSE
  )

  return(seurat_obj)
}


#' Full DoubletFinder Pipeline with Clustering and Singlet Filtering
#'
#' Comprehensive doublet detection workflow that includes normalization,
#' PCA, neighbour graph construction, clustering, and DoubletFinder
#' parameter sweep. Optionally returns only singlets.
#'
#' @param obj                  Seurat object.
#' @param etiqueta             Sample label for console messages.
#' @param PCs                  Integer vector of PCs (default 1:20).
#' @param resolution           Leiden resolution for pre-clustering (default 0.5).
#' @param return_singlets      If TRUE, subset to singlets before returning.
#' @param sct                  Whether the object was processed with SCTransform.
#' @return Seurat object with doublet_class metadata column.
#' @export
doubletfinder_pipeline <- function(obj,
                                   etiqueta       = "Sample",
                                   PCs            = 1:20,
                                   resolution     = 0.5,
                                   return_singlets = TRUE,
                                   sct            = FALSE) {
  message("Processing: ", etiqueta)

  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj)
  obj <- ScaleData(obj)
  obj <- RunPCA(obj, npcs = max(PCs))
  obj <- FindNeighbors(obj, dims = PCs)
  obj <- FindClusters(obj, resolution = resolution)

  sweep.res          <- paramSweep(obj, PCs = PCs, sct = sct)
  sweep.stats        <- summarizeSweep(sweep.res, GT = FALSE)
  sweep.stats$pK     <- as.numeric(as.character(sweep.stats$pK))
  sweep.stats$pN     <- as.numeric(as.character(sweep.stats$pN))

  best_row  <- sweep.stats[which.max(sweep.stats$BCreal), ]
  best.pK   <- best_row$pK
  best.pN   <- best_row$pN
  nExp      <- round(best.pN * ncol(obj))

  message("Best pK: ", best.pK, "  pN: ", best.pN, "  nExp: ", nExp)

  obj <- doubletFinder(
    obj,
    PCs        = PCs,
    pN         = best.pN,
    pK         = best.pK,
    nExp       = nExp,
    reuse.pANN = NULL,
    sct        = sct
  )

  # Flatten any data.frame metadata columns introduced by DoubletFinder
  for (col in colnames(obj@meta.data)) {
    if (is.data.frame(obj@meta.data[[col]])) {
      message("Fixing column: ", col)
      obj@meta.data[[col]] <- obj@meta.data[[col]][, 1]
    }
  }

  df_col           <- grep("DF.classifications", colnames(obj@meta.data), value = TRUE)
  obj$doublet_class <- obj[[df_col]]

  message("Summary for ", etiqueta, ":")
  print(table(obj$doublet_class))

  if (return_singlets) {
    obj <- subset(obj, subset = doublet_class == "Singlet")
    message("Retained singlets: ", ncol(obj))
  }

  return(obj)
}


#' Full Per-Sample Processing Pipeline
#'
#' Loads a CellBender HDF5 file, computes QC metrics, applies filters,
#' runs doublet detection, renames cells, and adds a condition column.
#'
#' Organellar gene patterns by organism:
#'   Arabidopsis thaliana : mt_pattern = "^ATMG",  cp_pattern = "^ATCG"
#'   Homo sapiens         : mt_pattern = "^MT-",   cp_pattern = NULL
#'   Mus musculus         : mt_pattern = "^mt-",   cp_pattern = NULL
#'   Zebrafish (D. rerio) : mt_pattern = "^mt-",   cp_pattern = NULL
#'
#' Set cp_pattern = NULL to skip chloroplast/plastid QC (non-plant organisms).
#'
#' @param sample_info  Named list with elements: file, label, condition.
#' @param min_features Minimum number of detected genes per cell.
#' @param max_features Maximum number of detected genes per cell.
#' @param min_counts   Minimum total UMI count per cell.
#' @param max_counts   Maximum total UMI count per cell.
#' @param max_mt       Maximum percentage of mitochondrial counts.
#' @param max_cp       Maximum percentage of chloroplast/plastid counts.
#'                     Ignored when cp_pattern = NULL.
#' @param mt_pattern   Regex pattern matching mitochondrial gene names.
#'                     Default "^ATMG" (Arabidopsis).
#' @param cp_pattern   Regex pattern matching chloroplast/plastid gene names,
#'                     or NULL to skip this metric. Default "^ATCG" (Arabidopsis).
#' @return Filtered Seurat object.
#' @export
process_sample <- function(sample_info,
                           min_features = 0,
                           max_features = Inf,
                           min_counts   = 0,
                           max_counts   = Inf,
                           max_mt       = 100,
                           max_cp       = 100,
                           mt_pattern   = "^ATMG",
                           cp_pattern   = "^ATCG") {

  obj <- load_cellbender_filtered_h5(sample_info$file, sample_info$label)

  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = mt_pattern)

  if (!is.null(cp_pattern)) {
    obj[["percent.cp"]] <- PercentageFeatureSet(obj, pattern = cp_pattern)
  }

  filter_expr <- quote(
    nFeature_RNA > min_features &
    nFeature_RNA < max_features &
    nCount_RNA   > min_counts   &
    nCount_RNA   < max_counts   &
    percent.mt   < max_mt
  )

  if (!is.null(cp_pattern)) {
    filter_expr <- bquote(.(filter_expr) & percent.cp < max_cp)
  }

  obj <- subset(obj, subset = eval(filter_expr))

  obj <- doubletfinder_pipeline(obj, etiqueta = sample_info$label)

  obj <- RenameCells(obj, add.cell.id = sample_info$condition)
  obj$condition <- sample_info$condition

  return(obj)
}


# =============================================================================
# 3. PSEUDOBULK
# =============================================================================

#' Generate Pseudobulk Count Matrices
#'
#' Aggregates single-cell counts by grouping variable (e.g., orig.ident) to
#' produce pseudobulk profiles. Optionally merges biological replicates by
#' condition using the $condition metadata column.
#'
#' @param seurat_obj       Seurat object with a $condition column.
#' @param group_by         Metadata column to aggregate by (default "orig.ident").
#' @param merge_replicates If TRUE, sum replicates sharing the same condition.
#' @return A matrix (group_by only) or a named list with:
#'   \describe{
#'     \item{by_sample}{Matrix with one column per group_by value.}
#'     \item{by_condition}{Matrix with one column per unique condition.}
#'   }
#' @export
generate_pseudobulk <- function(seurat_obj,
                                 group_by         = "orig.ident",
                                 merge_replicates = TRUE) {

  groups <- unique(seurat_obj@meta.data[[group_by]])
  cat("Generating pseudobulk for:", paste(groups, collapse = ", "), "\n")

  process_group <- function(group_name) {
    cells  <- subset(seurat_obj,
                     cells = colnames(seurat_obj)[
                       seurat_obj@meta.data[[group_by]] == group_name
                     ])
    layers <- grep("^counts", Layers(cells[["RNA"]]), value = TRUE)

    if (length(layers) == 0) {
      counts <- GetAssayData(cells[["RNA"]], layer = "data")
    } else if (length(layers) == 1) {
      counts <- GetAssayData(cells[["RNA"]], layer = layers)
    } else {
      mats   <- lapply(layers, function(x) GetAssayData(cells[["RNA"]], layer = x))
      counts <- Reduce(RowMergeSparseMatrices, mats)
    }

    gene_sums <- Matrix::rowSums(counts)
    cat(" ", group_name, "->", ncol(cells), "cells,", length(gene_sums), "genes\n")
    return(gene_sums)
  }

  pseudobulk_list <- lapply(groups, process_group)
  names(pseudobulk_list) <- groups

  all_genes         <- unique(unlist(lapply(pseudobulk_list, names)))
  pseudobulk_matrix <- sapply(pseudobulk_list, function(x) {
    v       <- x[all_genes]
    v[is.na(v)] <- 0
    return(v)
  })
  rownames(pseudobulk_matrix) <- all_genes

  if (merge_replicates) {
    conditions <- unique(seurat_obj$condition)

    merged_matrix <- sapply(conditions, function(cond) {
      cols <- grep(paste0("^", cond), colnames(pseudobulk_matrix), value = TRUE)
      if (length(cols) == 0) {
        cols <- colnames(pseudobulk_matrix)[grepl(cond, colnames(pseudobulk_matrix))]
      }
      if (length(cols) == 1) return(pseudobulk_matrix[, cols])
      return(rowSums(pseudobulk_matrix[, cols, drop = FALSE]))
    })
    colnames(merged_matrix) <- conditions

    cat("\nMerged replicates:\n")
    print(colnames(merged_matrix))

    return(list(
      by_sample    = pseudobulk_matrix,
      by_condition = merged_matrix
    ))
  }

  return(pseudobulk_matrix)
}


#' Normalize Pseudobulk and Bulk Counts with DESeq2
#'
#' Applies DESeq2 size-factor normalization to a pair of pseudobulk and
#' bulk count vectors sharing common genes, then log2-transforms the result.
#'
#' @param pseudobulk_counts Named numeric vector of pseudobulk raw counts.
#' @param bulk_counts        Named numeric vector of bulk raw counts.
#' @return Data frame with columns: gene, pseudobulk, bulk (log2-normalized).
#' @export
normalizar_bulk_pseudobulk <- function(pseudobulk_counts, bulk_counts) {

  common_genes <- intersect(names(pseudobulk_counts), names(bulk_counts))
  if (length(common_genes) < 10) stop("Too few common genes.")

  counts_matrix <- data.frame(
    pseudobulk = round(pseudobulk_counts[common_genes]),
    bulk       = round(bulk_counts[common_genes])
  )
  rownames(counts_matrix) <- common_genes

  condition <- factor(c("pseudobulk", "bulk"))
  col_data  <- data.frame(condition = condition)

  dds        <- DESeqDataSetFromMatrix(countData = counts_matrix,
                                       colData   = col_data,
                                       design    = ~ condition)
  dds        <- estimateSizeFactors(dds)
  norm_counts <- counts(dds, normalized = TRUE)

  log_norm_counts <- log2(norm_counts + 1)

  df <- data.frame(
    gene       = rownames(log_norm_counts),
    pseudobulk = log_norm_counts[, "pseudobulk"],
    bulk       = log_norm_counts[, "bulk"]
  )

  return(df)
}


#' Classify Genes by Linear Model Residuals
#'
#' Fits a linear model (bulk ~ pseudobulk) and classifies each gene as
#' Upregulated, Downregulated, or Consistent based on residual magnitude.
#'
#' @param df     Data frame with columns pseudobulk and bulk.
#' @param umbral Absolute residual threshold for classification (default 5).
#' @return The input data frame with additional columns: residuals, status.
#' @export
clasificar_residuos <- function(df, umbral = 5) {

  modelo      <- lm(bulk ~ pseudobulk, data = df)
  df$residuals <- resid(modelo)

  df$status <- case_when(
    df$residuals >  umbral ~ "Upregulated",
    df$residuals < -umbral ~ "Downregulated",
    TRUE                   ~ "Consistent"
  )

  return(df)
}


# =============================================================================
# 4. CLUSTERING UTILITIES
# =============================================================================

#' Unify Cluster Ident Names
#'
#' Removes numeric suffixes (e.g., ".1", "_2") appended to cluster names
#' during integration or merging, restoring the canonical label.
#'
#' @param obj Seurat object.
#' @return Seurat object with updated Idents.
#' @export
unificar_nombres <- function(obj) {
  old_levels <- levels(obj)
  new_levels <- gsub("[._][0-9]+$", "", old_levels)
  new_ids    <- setNames(new_levels, old_levels)
  obj        <- RenameIdents(obj, new_ids)
  return(obj)
}


#' Merge Multiple Count Layers into One Sparse Matrix
#'
#' Combines multiple sparse count layers from the RNA assay of a Seurat 5
#' object (e.g., after merging samples) into a single dgCMatrix.
#'
#' @param obj   Seurat object.
#' @param capas Character vector of layer names to merge.
#' @return A sparse dgCMatrix of merged counts.
#' @export
unir_layers_counts <- function(obj, capas) {

  if (length(capas) == 1) {
    return(GetAssayData(obj[["RNA"]], layer = capas))
  }

  message("Merging ", length(capas), " layers...")
  mats   <- lapply(capas, function(x) GetAssayData(obj[["RNA"]], layer = x))
  merged <- Reduce(RowMergeSparseMatrices, mats)
  return(merged)
}


# =============================================================================
# 5. MARKER-BASED AND REFERENCE-BASED ANNOTATION
# =============================================================================

#' Find and Cache All Marker Genes
#'
#' Runs FindAllMarkers on the seurat_clusters identity and writes results to a
#' TSV file. If the output file already exists (and force = FALSE), it loads
#' the cached result instead.
#'
#' @param seurat_obj      Seurat object.
#' @param output_file     Path to output TSV file.
#' @param only_pos        Return only positive markers (default TRUE).
#' @param min_pct         Minimum fraction of cells expressing the gene.
#' @param logfc_threshold Minimum log2 fold change threshold.
#' @param force           If TRUE, recompute even if output_file exists.
#' @return Data frame of marker genes.
#' @export
find_markers <- function(seurat_obj,
                         output_file     = "results/FindAllMarkers.tsv",
                         only_pos        = TRUE,
                         min_pct         = 0.25,
                         logfc_threshold = 0.25,
                         force           = FALSE) {

  seurat_obj      <- JoinLayers(seurat_obj)
  Idents(seurat_obj) <- "seurat_clusters"

  if (file.exists(output_file) && !force) {
    cat("Loading cached markers:", output_file, "\n")
    markers <- read.table(output_file, header = TRUE, sep = "\t", quote = "")
  } else {
    cat("Computing markers...\n")
    markers <- FindAllMarkers(
      seurat_obj,
      only.pos        = only_pos,
      min.pct         = min_pct,
      logfc.threshold = logfc_threshold
    )
    write.table(markers, output_file, quote = FALSE, sep = "\t", row.names = FALSE)
    cat("Saved to:", output_file, "\n")
  }

  return(markers)
}


#' Annotate Clusters Using a Known Marker Gene List
#'
#' Crosses FindAllMarkers output against a user-provided reference table
#' (columns: gene, cell.types) to assign a cell-type label to each cluster.
#' The best-matching cell type per cluster is selected by adjusted p-value.
#'
#' @param seurat_obj     Seurat object.
#' @param markers        Data frame returned by find_markers().
#' @param reference_file Path to a tab-delimited file with columns gene and
#'                       cell.types. If NULL, a file-chooser dialog opens.
#' @return Seurat object with a celltype metadata column.
#' @export
annotate_by_markers <- function(seurat_obj,
                                markers,
                                reference_file = NULL) {

  if (is.null(reference_file)) {
    reference_file <- file.choose(caption = "Select reference file (gene | cell.types)")
  }
  cat("Using reference:", reference_file, "\n")

  reference <- read.table(reference_file, header = TRUE, sep = "\t", quote = "")

  merged <- merge(markers, reference, by.x = "gene", by.y = "gene")
  merged <- merged[order(merged$cluster, merged$p_val_adj), ]
  merged <- merged[!duplicated(merged$cluster), ]

  cat("\nMatches found:\n")
  print(merged[, c("cluster", "gene", "cell.types")])

  Idents(seurat_obj) <- "seurat_clusters"
  new_ids             <- merged$cell.types
  names(new_ids)      <- merged$cluster
  seurat_obj          <- RenameIdents(seurat_obj, new_ids)
  seurat_obj$celltype <- Idents(seurat_obj)

  cat("\nFinal annotation:\n")
  print(table(seurat_obj$celltype))

  return(seurat_obj)
}


#' Annotate Clusters via Label Transfer from a Reference Seurat Object
#'
#' Uses Seurat's FindTransferAnchors / TransferData workflow to project
#' cell-type labels from a reference dataset onto the query object.
#'
#' @param seurat_obj   Query Seurat object.
#' @param reference_obj Reference Seurat object. If NULL, a file-chooser
#'                     dialog prompts for an RDS file.
#' @param reference_col Metadata column in reference_obj containing labels.
#'                     If NULL, the user selects interactively.
#' @param dims         Integer vector of PCs for anchor finding (default 1:30).
#' @return Query Seurat object with a celltype_reference metadata column.
#' @export
annotate_by_reference <- function(seurat_obj,
                                  reference_obj = NULL,
                                  reference_col = NULL,
                                  dims          = 1:30) {

  if (is.null(reference_obj)) {
    ref_file      <- file.choose(caption = "Select reference Seurat object (.rds)")
    cat("Loading reference:", ref_file, "\n")
    reference_obj <- readRDS(ref_file)
  }

  if (is.null(reference_col)) {
    cat("\nAvailable columns in reference:\n")
    cols <- colnames(reference_obj@meta.data)
    for (i in seq_along(cols)) cat(" ", i, "->", cols[i], "\n")
    selection     <- as.integer(readline("Select column number: "))
    reference_col <- cols[selection]
  }

  cat("Using column:", reference_col, "\n")

  anchors <- FindTransferAnchors(
    reference = reference_obj,
    query     = seurat_obj,
    dims      = dims
  )

  predictions <- TransferData(
    anchorset = anchors,
    refdata   = reference_obj@meta.data[[reference_col]],
    dims      = dims
  )

  seurat_obj$celltype_reference <- predictions$predicted.id

  cat("\nReference-based annotation:\n")
  print(table(seurat_obj$celltype_reference))

  return(seurat_obj)
}


# =============================================================================
# 6. EXPORT
# =============================================================================

#' Export Seurat Object to Scanpy-Compatible h5ad Format
#'
#' Converts a Seurat object to SingleCellExperiment and writes an AnnData
#' h5ad file. Tries zellkonverter first; falls back to SeuratDisk if needed.
#' Counts and logcounts assays are included, together with selected
#' dimensionality reduction embeddings.
#'
#' @param seurat_obj  Seurat object to export.
#' @param outfile     Output file path (should end in .h5ad).
#' @param assay_name  Assay to export (default "RNA").
#' @param use_reduc   Character vector of reductions to include
#'                    (default c("pca", "umap", "harmony")).
#' @param X_name      Name of the .X slot in the AnnData object (default
#'                    "logcounts").
#' @param overwrite   If TRUE, overwrite an existing file (default TRUE).
#' @return Invisibly, the SingleCellExperiment intermediate.
#' @export
exportar_para_scanpy <- function(seurat_obj,
                                 outfile,
                                 assay_name = "RNA",
                                 use_reduc  = c("pca", "umap", "harmony"),
                                 X_name     = "logcounts",
                                 overwrite  = TRUE) {

  stopifnot(inherits(seurat_obj, "Seurat"))
  if (!dir.exists(dirname(outfile))) dir.create(dirname(outfile), recursive = TRUE)

  message("Converting Seurat -> SCE...")
  sce <- as.SingleCellExperiment(seurat_obj, assay = assay_name)

  if (is.null(rownames(sce))) stop("Object missing rownames.")
  rownames(sce) <- make.unique(rownames(sce))

  if (is.null(colnames(sce)) || any(!nzchar(colnames(sce)))) {
    colnames(sce) <- paste0("cell", seq_len(ncol(sce)))
  }
  stopifnot(!anyDuplicated(colnames(sce)))

  message("Extracting counts and logcounts...")
  assay(sce, "counts")    <- Seurat::GetAssayData(seurat_obj, assay = assay_name,
                                                   slot = "counts")
  assay(sce, "logcounts") <- Seurat::GetAssayData(seurat_obj, assay = assay_name,
                                                   slot = "data")

  if (ncol(colData(sce)) == 0) {
    colData(sce) <- S4Vectors::DataFrame(row.names = colnames(sce))
  }
  colData(sce) <- cbind(
    colData(sce),
    seurat_obj@meta.data[colnames(sce), , drop = FALSE]
  )

  message("Exporting reductions...")
  reds <- Seurat::Reductions(seurat_obj)
  for (red in use_reduc) {
    if (red %in% reds) {
      message("  Including: ", red)
      reducedDims(sce)[[toupper(red)]] <- seurat_obj@reductions[[red]]@cell.embeddings
    }
  }

  if (file.exists(outfile) && overwrite) file.remove(outfile)

  ok <- FALSE
  if (requireNamespace("zellkonverter", quietly = TRUE)) {
    message("Writing h5ad with zellkonverter...")
    zellkonverter::writeH5AD(sce, file = outfile, X_name = X_name)
    ok <- TRUE
  } else if (requireNamespace("SeuratDisk", quietly = TRUE)) {
    message("Falling back to SeuratDisk...")
    tmp_h5seu <- file.path(tempdir(), paste0(basename(outfile), ".h5seurat"))
    if (file.exists(tmp_h5seu)) file.remove(tmp_h5seu)
    SeuratDisk::SaveH5Seurat(seurat_obj, filename = tmp_h5seu, overwrite = TRUE)
    SeuratDisk::Convert(source = tmp_h5seu, dest = "h5ad", overwrite = TRUE)
    gen_h5ad <- sub("\\.h5seurat$", ".h5ad", tmp_h5seu)
    if (!file.exists(gen_h5ad)) stop("Failed to generate h5ad.")
    file.rename(gen_h5ad, outfile)
    ok <- TRUE
  } else {
    stop("Install 'zellkonverter' or 'SeuratDisk' to export h5ad.")
  }

  if (!ok || !file.exists(outfile)) stop("Export failed: ", outfile)
  message("Export complete: ", outfile)
  invisible(sce)
}


# =============================================================================
# END OF FUNCTIONS
# =============================================================================
