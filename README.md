# ScRNASeq-Docker — Helper Function Library

This repository contains the pipeline helper scripts used in the **Methods in Molecular Biology** protocol for single-cell RNA-seq analysis of *Arabidopsis thaliana* (readily extensible to other organisms). The library covers the complete analytical trajectory: CellBender-aware data loading, QC visualization, doublet detection, pseudobulk aggregation, differential expression with DESeq2, and GO enrichment — all orchestrated through a set of modular, well-documented R functions.

---

## Repository files

| File | Purpose |
|---|---|
| `load_libraries.R` | Loads and, if necessary, installs all R package dependencies |
| `ScRNA_Analysis_Functions.R` | Core helper function library (documented in full below) |
| `custom_seurat.R` | Project-specific Seurat extensions and theme overrides |
| `scrnaseq_pipeline.R` | Top-level pipeline script that calls the helper functions end-to-end |
| `analysis.Rmd` / `analysis_code.pdf` | Reproducible analysis notebook and its rendered PDF for the protocol |

---

## Function reference

### 1. QC and Visualization

---

#### `load_cellbender_filtered_h5()`

Reads a filtered expression matrix from a CellBender HDF5 output file and returns a Seurat object containing the raw counts.

| Parameter | Type | Description |
|---|---|---|
| `h5_path` | `character` | Path to the filtered `.h5` file produced by CellBender |
| `project` | `character` | Project name stored in Seurat object metadata (default `"Sample"`) |

```r
seu <- load_cellbender_filtered_h5("results/sample_filtered.h5", project = "Root_WT")
```

---

#### `plot_qc_violin_grid()`

Produces a violin plot grid showing `nFeature_RNA`, `nCount_RNA`, `percent.mt`, and (if present) `percent.cp` for a single Seurat object, with the cell count in the plot title.

| Parameter | Type | Description |
|---|---|---|
| `obj1` | `Seurat` | Seurat object to visualize |
| `label` | `character` | Condition label displayed as the group and plot title |
| `color` | `character` | Fill color for the violin plots |

```r
p <- plot_qc_violin_grid(seu, label = "WT_rep1", color = "#66c2a5")
```

---

#### `resumen_nFeature_plot()`

Creates a boxplot of `nFeature_RNA` distributions across a list of Seurat objects, alongside printed quartile and quintile summary tables rendered as grid graphics.

| Parameter | Type | Description |
|---|---|---|
| `obj_list` | `list` | List of Seurat objects |
| `etiquetas` | `character` | Labels for each object (default: `"Group1"`, `"Group2"`, ...) |
| `colores` | `character` | Named or positional color vector; auto-assigned from ColorBrewer if `NULL` |

```r
resumen_nFeature_plot(list(wt, mut), etiquetas = c("WT", "Mutant"), colores = c("steelblue", "tomato"))
```

---

### 2. Preprocessing and Doublet Detection

---

#### `preprocesar_y_doubletfinder()`

Normalizes a Seurat object, finds variable features, scales, runs PCA, performs the DoubletFinder parameter sweep, and returns the object with doublet classifications added to metadata.

| Parameter | Type | Description |
|---|---|---|
| `seurat_obj` | `Seurat` | Input Seurat object |
| `pcs` | `integer` | PC range to use (e.g. `1:20`) |
| `expected_doublet_rate` | `numeric` | Expected fraction of doublets (default `0.075`) |
| `project_id` | `character` | Label for log messages (default `"sample"`) |

```r
seu <- preprocesar_y_doubletfinder(seu, pcs = 1:20, expected_doublet_rate = 0.08, project_id = "Root_WT")
```

---

#### `doubletfinder_pipeline()`

Comprehensive doublet detection pipeline that adds neighborhood graph construction and clustering before the DoubletFinder parameter sweep, with optional automatic filtering to singlets only.

| Parameter | Type | Description |
|---|---|---|
| `obj` | `Seurat` | Input Seurat object |
| `etiqueta` | `character` | Sample label for progress messages (default `"Sample"`) |
| `PCs` | `integer` | PC range for all dimension-reduction steps (default `1:20`) |
| `resolution` | `numeric` | Louvain/Leiden clustering resolution (default `0.5`) |
| `return_singlets` | `logical` | If `TRUE`, subset to singlets before returning (default `TRUE`) |
| `sct` | `logical` | Whether SCTransform normalization was used (default `FALSE`) |

```r
seu_clean <- doubletfinder_pipeline(seu, etiqueta = "Root_WT", PCs = 1:30, return_singlets = TRUE)
```

---

#### `load_sample()`

Loads a CellBender HDF5 file, computes mitochondrial and (optionally) chloroplast percentages, prefixes cell barcodes with the condition label, and returns the annotated object without applying any filters — useful for inspecting raw QC metrics before deciding thresholds.

| Parameter | Type | Description |
|---|---|---|
| `sample_info` | `list` | Named list with fields `file` (path), `label` (project name), and `condition` (barcode prefix) |
| `mt_pattern` | `character` | Regex matching mitochondrial gene IDs (default `"^ATMG"`) |
| `cp_pattern` | `character` | Regex matching chloroplast gene IDs; `NULL` to skip (default `"^ATCG"`) |

```r
seu_raw <- load_sample(list(file = "sample.h5", label = "WT", condition = "WT"), mt_pattern = "^ATMG", cp_pattern = "^ATCG")
```

---

#### `filter_sample()`

Applies QC thresholds to an already-annotated Seurat object (output of `load_sample()`) and optionally runs the full DoubletFinder pipeline on the filtered cells.

| Parameter | Type | Description |
|---|---|---|
| `obj` | `Seurat` | Annotated Seurat object (from `load_sample()`) |
| `min_features` | `numeric` | Minimum `nFeature_RNA` (default `200`) |
| `max_features` | `numeric` | Maximum `nFeature_RNA` (default `Inf`) |
| `min_counts` | `numeric` | Minimum `nCount_RNA` (default `0`) |
| `max_counts` | `numeric` | Maximum `nCount_RNA` (default `Inf`) |
| `max_mt` | `numeric` | Maximum mitochondrial percent (default `5`) |
| `max_cp` | `numeric` | Maximum chloroplast percent; ignored if `percent.cp` is absent (default `100`) |
| `run_doubletfinder` | `logical` | Whether to run DoubletFinder after filtering (default `TRUE`) |

```r
seu_filt <- filter_sample(seu_raw, min_features = 300, max_features = 8000, max_mt = 5)
```

---

#### `process_sample()`

Convenience wrapper that calls `load_sample()` followed by `filter_sample()` in a single step; use this when you do not need to inspect raw QC plots before setting thresholds.

| Parameter | Type | Description |
|---|---|---|
| `sample_info` | `list` | Named list with fields `file`, `label`, and `condition` (see `load_sample()`) |
| `mt_pattern` | `character` | Mitochondrial gene regex (default `"^ATMG"`) |
| `cp_pattern` | `character` | Chloroplast gene regex; `NULL` to skip (default `"^ATCG"`) |
| `min_features` | `numeric` | Minimum `nFeature_RNA` (default `200`) |
| `max_features` | `numeric` | Maximum `nFeature_RNA` (default `Inf`) |
| `min_counts` | `numeric` | Minimum `nCount_RNA` (default `0`) |
| `max_counts` | `numeric` | Maximum `nCount_RNA` (default `Inf`) |
| `max_mt` | `numeric` | Maximum mitochondrial percent (default `5`) |
| `max_cp` | `numeric` | Maximum chloroplast percent (default `100`) |
| `run_doubletfinder` | `logical` | Whether to run DoubletFinder (default `TRUE`) |

```r
seu <- process_sample(list(file = "sample.h5", label = "WT", condition = "WT"), min_features = 300, max_mt = 5)
```

---

### 3. Bulk / Pseudobulk Utilities

---

#### `normalizar_bulk_pseudobulk()`

Applies DESeq2 size-factor normalization followed by log2 transformation to a pair of pseudobulk and bulk count vectors, restricting the analysis to their common gene set.

| Parameter | Type | Description |
|---|---|---|
| `pseudobulk_counts` | `named numeric` | Named vector of pseudobulk raw counts (genes as names) |
| `bulk_counts` | `named numeric` | Named vector of bulk RNA-seq raw counts (genes as names) |

```r
df_norm <- normalizar_bulk_pseudobulk(pseudo_vec, bulk_vec)
```

---

#### `clasificar_residuos()`

Fits a linear model of bulk expression on pseudobulk expression and classifies each gene as `"Upregulated"`, `"Downregulated"`, or `"Consistent"` based on whether its residual exceeds the specified threshold.

| Parameter | Type | Description |
|---|---|---|
| `df` | `data.frame` | Data frame with columns `pseudobulk` and `bulk` (e.g. from `normalizar_bulk_pseudobulk()`) |
| `umbral` | `numeric` | Absolute residual threshold for classification (default `5`) |

```r
df_classified <- clasificar_residuos(df_norm, umbral = 5)
```

---

#### `generate_pseudobulk()`

Aggregates single-cell counts by a metadata grouping variable to produce a pseudobulk count matrix; optionally merges per-sample columns into per-condition totals.

| Parameter | Type | Description |
|---|---|---|
| `seurat_obj` | `Seurat` | Seurat object containing raw counts |
| `group_by` | `character` | Metadata column name used to group cells (default `"orig.ident"`) |
| `merge_replicates` | `logical` | If `TRUE`, return both per-sample and per-condition matrices (default `TRUE`) |

```r
pb <- generate_pseudobulk(seu, group_by = "orig.ident", merge_replicates = TRUE)
# Access: pb$by_sample, pb$by_condition
```

---

#### `plot_replicate_correlation()`

Computes pairwise Pearson correlations across columns of a pseudobulk matrix and renders a `pheatmap` of the correlation matrix with numeric annotations.

| Parameter | Type | Description |
|---|---|---|
| `pseudobulk_mat` | `matrix` | Numeric genes-by-samples matrix (e.g. from `generate_pseudobulk()` or `hacer_pseudobulk()`) |
| `main` | `character` | Heatmap title (default `"Replicate Correlation"`) |

```r
plot_replicate_correlation(pb$by_sample, main = "QC: Replicate Correlation")
```

---

### 4. Seurat Utilities

---

#### `unificar_nombres()`

Removes numeric suffixes (e.g. `.1`, `_2`) from cluster identity level names, unifying duplicated cluster labels that can arise after merging Seurat objects.

| Parameter | Type | Description |
|---|---|---|
| `obj` | `Seurat` | Seurat object whose active ident levels should be cleaned |

```r
seu <- unificar_nombres(seu)
```

---

#### `mostrar_tabla()`

Creates and renders a side-by-side cell-type count comparison table (filtered vs. CellBender annotations) using grid graphics.

| Parameter | Type | Description |
|---|---|---|
| `filtered_vec` | `character` | Annotation vector from the filtered object |
| `cellbender_vec` | `character` | Annotation vector from the CellBender object |
| `titulo` | `character` | Table title (default `"Annotations"`) |

```r
mostrar_tabla(seu_filt$celltype, seu_raw$celltype, titulo = "Cell-type counts")
```

---

#### `exportar_para_scanpy()`

Converts a Seurat object to SingleCellExperiment and writes it as an `.h5ad` file compatible with Scanpy/AnnData; prefers `zellkonverter` but falls back to `SeuratDisk` if necessary.

| Parameter | Type | Description |
|---|---|---|
| `seurat_obj` | `Seurat` | Seurat object to export |
| `outfile` | `character` | Output file path (must end in `.h5ad`) |
| `assay_name` | `character` | Assay to export (default `"RNA"`) |
| `use_reduc` | `character` | Reductions to embed in the h5ad (default `c("pca","umap","harmony")`) |
| `X_name` | `character` | Assay layer stored as `.X` in Scanpy (default `"logcounts"`) |
| `overwrite` | `logical` | Overwrite an existing file (default `TRUE`) |

```r
exportar_para_scanpy(seu, outfile = "results/atlas.h5ad", use_reduc = c("pca", "umap"))
```

---

#### `safe_vln()`

Thin wrapper around `VlnPlot` that groups cells by `orig.ident`, suppresses individual points, and applies a custom fill palette — safe for use inside RMarkdown without interactive prompts.

| Parameter | Type | Description |
|---|---|---|
| `obj` | `Seurat` | Seurat object |
| `feature` | `character` | Gene name or metadata column to plot |
| `colors` | `character` | Named or positional color palette |

```r
p <- safe_vln(seu, feature = "AT1G01060", colors = sample_colors)
```

---

#### `unir_layers_counts()`

Merges multiple sparse count matrices from different RNA assay layers into a single sparse matrix, handling the single-layer case without unnecessary copying.

| Parameter | Type | Description |
|---|---|---|
| `obj` | `Seurat` | Seurat object containing the RNA assay |
| `capas` | `character` | Character vector of layer names to merge |

```r
merged_mat <- unir_layers_counts(seu, capas = c("counts.WT", "counts.Mut"))
```

---

### 5. Annotation

---

#### `find_markers()`

Runs `FindAllMarkers` on Seurat clusters and writes results to a TSV cache file, loading from the cache on subsequent calls unless `force = TRUE`.

| Parameter | Type | Description |
|---|---|---|
| `seurat_obj` | `Seurat` | Seurat object with `seurat_clusters` active identity |
| `output_file` | `character` | Path for the TSV cache (default `"results/FindAllMarkers.tsv"`) |
| `only_pos` | `logical` | Return only positive markers (default `TRUE`) |
| `min_pct` | `numeric` | Minimum fraction of cells expressing the gene (default `0.25`) |
| `logfc_threshold` | `numeric` | Log fold-change threshold (default `0.25`) |
| `force` | `logical` | Recompute even if a cache file exists (default `FALSE`) |

```r
markers <- find_markers(seu, output_file = "results/FindAllMarkers.tsv", force = FALSE)
```

---

#### `annotate_by_markers()`

Crosses `FindAllMarkers` output with a two-column reference table (`gene | cell.types`) to assign the best-matching cell-type label to each cluster, then stores the result in the `celltype` metadata column.

| Parameter | Type | Description |
|---|---|---|
| `seurat_obj` | `Seurat` | Seurat object |
| `markers` | `data.frame` | Output of `find_markers()` |
| `reference_file` | `character` | Path to the tab-separated reference table; a file chooser dialog is shown if `NULL` |

```r
seu <- annotate_by_markers(seu, markers, reference_file = "refs/cell_type_markers.tsv")
```

---

#### `annotate_by_reference()`

Transfers cell-type labels from a reference Seurat object to the query using `FindTransferAnchors` and `TransferData`, storing predictions in a `celltype_reference` metadata column.

| Parameter | Type | Description |
|---|---|---|
| `seurat_obj` | `Seurat` | Query Seurat object |
| `reference_obj` | `Seurat` | Reference Seurat object; a file chooser is shown if `NULL` |
| `reference_col` | `character` | Metadata column in the reference holding cell-type labels; interactive selection if `NULL` |
| `dims` | `integer` | Dimensions for anchor finding (default `1:30`) |

```r
seu <- annotate_by_reference(seu, reference_obj = ref_atlas, reference_col = "celltype", dims = 1:30)
```

---

#### `subclustar_tipo()`

Subsets a Seurat object to one or more named cell types, then re-runs PCA, UMAP, neighbor finding, and clustering at the specified resolution to reveal sub-populations within that compartment.

| Parameter | Type | Description |
|---|---|---|
| `obj` | `Seurat` | Seurat object with cell-type annotations |
| `tipo` | `character` | Cell type(s) to subset (must match values in `annot_col`) |
| `annot_col` | `character` | Metadata column holding cell-type labels (default `"annotation_agrupada"`) |
| `resolution` | `numeric` | Clustering resolution (default `0.3`) |
| `dims` | `integer` | Dimensions for UMAP and neighbor finding (default `1:20`) |

```r
guard_sub <- subclustar_tipo(seu, tipo = "Guard Cell", annot_col = "celltype", resolution = 0.3)
```

---

### 6. Pseudobulk, DESeq2, Volcano, Heatmap

---

#### `asignar_pseudoreplicados()`

Randomly assigns cells within each condition to a specified number of pseudo-replicate groups, enabling DESeq2-based differential expression on experiments without biological replicates.

| Parameter | Type | Description |
|---|---|---|
| `obj` | `Seurat` | Seurat object with an `orig.ident_uni` metadata column |
| `condiciones` | `character` | Conditions to include; `NULL` uses all detected conditions |
| `n_reps` | `integer` | Number of pseudo-replicates per condition (default `3`) |
| `seed` | `integer` | Random seed for reproducibility (default `1807`) |

```r
seu <- asignar_pseudoreplicados(seu, condiciones = c("WT", "Mutant"), n_reps = 3, seed = 42)
```

---

#### `hacer_pseudobulk()`

Aggregates counts per pseudo-replicate group using `AggregateExpression` and returns a tidy data frame (genes as rows, replicates as columns) ready for DESeq2 input.

| Parameter | Type | Description |
|---|---|---|
| `obj` | `Seurat` | Seurat object with a `replicate` metadata column (from `asignar_pseudoreplicados()`) |

```r
counts_mat <- hacer_pseudobulk(seu)
```

---

#### `correr_deseq2()`

Builds a DESeqDataSet from a pseudobulk count matrix, auto-detects condition levels from column names, runs DESeq2, and writes per-comparison result CSV files to disk.

| Parameter | Type | Description |
|---|---|---|
| `counts_mat` | `matrix` | Integer genes-by-samples count matrix |
| `comparaciones` | `list` | List of named lists, each with `conds` (length-2 character: reference then treatment) and `tag` (output label) |
| `output_dir` | `character` | Root directory; results written to `output_dir/tag/DESeq2_tag.csv` |
| `tipo` | `character` | Optional prefix added to output filenames (default `NULL`) |

```r
correr_deseq2(counts_mat,
              comparaciones = list(list(conds = c("WT", "Mutant"), tag = "WT_vs_Mutant")),
              output_dir = "results/DESeq2")
```

---

#### `hacer_volcano()`

Reads a DESeq2 CSV output file and produces a volcano plot colored by significance category (upregulated / downregulated / not significant) with user-defined fold-change and p-value cutoff lines.

| Parameter | Type | Description |
|---|---|---|
| `file` | `character` | Path to the DESeq2 results CSV file |
| `padj_cut` | `numeric` | Adjusted p-value significance cutoff (default `0.05`) |
| `lfc_cut` | `numeric` | Log2 fold-change cutoff (default `1`) |

```r
p <- hacer_volcano("results/DESeq2/WT_vs_Mutant/DESeq2_WT_vs_Mutant.csv", padj_cut = 0.05, lfc_cut = 1)
```

---

#### `procesar_deseq2_resultado()`

Reads a DESeq2 CSV, classifies each gene as up (`1`), down (`-1`), or unchanged (`0`), extracts log fold-change values for significant genes, and writes a filtered significant-gene CSV to disk.

| Parameter | Type | Description |
|---|---|---|
| `file_path` | `character` | Path to the DESeq2 results CSV file |
| `output_dir` | `character` | Directory where the filtered CSV is written |
| `padj_cut` | `numeric` | Adjusted p-value cutoff (default `0.05`) |
| `lfc_cut` | `numeric` | Log2 fold-change cutoff (default `1`) |

```r
res <- procesar_deseq2_resultado("results/DESeq2/WT_vs_Mutant/DESeq2_WT_vs_Mutant.csv",
                                  output_dir = "results/DESeq2/filtered")
# Returns: res$class and res$logfc data frames
```

---

#### `hacer_heatmap()`

Renders a hierarchically clustered heatmap (rows clustered by Euclidean distance with dynamic tree cut, columns by PCA-based distance) with per-gene cluster color annotations and a blue-black-yellow log2FC color scale.

| Parameter | Type | Description |
|---|---|---|
| `matriz` | `matrix` | Numeric genes-by-conditions matrix (e.g. log2FC values) |
| `min_genes` | `integer` | Minimum cluster size for dynamic tree cut (default `1`) |
| `deepSplit_val` | `integer` | `deepSplit` parameter for `cutreeDynamic` (default `0`) |
| `breaks` | `numeric` | Two-element `c(min, max)` vector for the color scale (default `c(-5, 5)`) |

```r
hacer_heatmap(lfc_matrix, min_genes = 5, deepSplit_val = 1, breaks = c(-3, 3))
```

---

#### `hacer_dotplot_marcadores()`

Builds a coordinate-flipped `DotPlot` where cell types and marker genes follow user-defined orders, producing a near-diagonal expression pattern useful for cell-type validation figures; optionally saves the plot as a PDF.

| Parameter | Type | Description |
|---|---|---|
| `seurat_obj` | `Seurat` | Seurat object with cell-type annotations |
| `marks` | `data.frame` | Data frame with columns `gene` and `cell.types` |
| `annot_col` | `character` | Metadata column holding cell-type labels (default `"celltype_reference_curated"`) |
| `cell_order` | `character` | Desired top-to-bottom cell-type order; unlisted types appended at the end |
| `clusters_remove` | `character` | Cell-type labels to exclude from the plot (default `NULL`) |
| `rename_map` | `named character` | Optional mapping to rename cell types before plotting |
| `outfile` | `character` | PDF output path; `NULL` skips saving (default `NULL`) |
| `width` | `numeric` | PDF width in inches (default `20`) |
| `height` | `numeric` | PDF height in inches (default `10`) |
| `dot_scale` | `numeric` | Dot size scaling factor (default `12`) |
| `base_size` | `numeric` | Base font size (default `18`) |

```r
p <- hacer_dotplot_marcadores(seu, marks = marker_df, annot_col = "celltype",
                               cell_order = c("Epidermis", "Vasculature", "Mesophyll"),
                               outfile = "figures/dotplot_markers.pdf")
```

---

### 7. GO Enrichment

---

#### `correr_enriquecimiento_go()`

Iterates over columns of a binary classification matrix, runs `clusterProfiler::enrichGO` for each set of upregulated genes, writes raw and gene-symbol-readable result tables to disk, and optionally simplifies redundant GO terms.

| Parameter | Type | Description |
|---|---|---|
| `tabla` | `matrix` | Binary genes-by-comparisons matrix; rows with value `1` are tested |
| `universo` | `character` | Background gene ID vector for enrichment testing |
| `espacio` | `character` | GO namespace: `"BP"`, `"MF"`, or `"CC"` |
| `orgdb` | `OrgDb` | OrgDb annotation object (default `org.At.tair.db`) |
| `keytype` | `character` | Key type matching rownames of `tabla` (default `"TAIR"`) |
| `qvalueCutoff` | `numeric` | Q-value significance cutoff (default `0.05`) |
| `pvalueCutoff` | `numeric` | P-value significance cutoff (default `0.05`) |
| `simplificar` | `logical` | If `TRUE`, remove redundant GO terms with `simplify()` (default `FALSE`) |
| `umbral_simply` | `numeric` | Similarity cutoff passed to `simplify()` (default `0.7`) |
| `output_dir` | `character` | Directory for output text files (default `"results/Enrichment"`) |

```r
go_results <- correr_enriquecimiento_go(class_matrix, universo = all_genes, espacio = "BP",
                                         orgdb = org.At.tair.db, keytype = "TAIR",
                                         output_dir = "results/GO")
```

---

#### `podar_go()`

Applies `gofilter` to each `enrichResult` in a named list to retain only GO terms at or below a specified ontology level, writing filtered tables to disk.

| Parameter | Type | Description |
|---|---|---|
| `resuGO` | `list` | Named list of `enrichResult` objects (from `correr_enriquecimiento_go()`) |
| `nivel` | `integer` | Maximum GO level to retain |
| `espacio` | `character` | GO namespace string (used in output filenames) |
| `qvalueCutoff` | `numeric` | Q-value cutoff (used in output filenames) |
| `simplificar` | `logical` | Affects the output filename suffix (default `FALSE`) |
| `output_dir` | `character` | Directory for output files (default `"results/Enrichment"`) |

```r
go_pruned <- podar_go(go_results, nivel = 4, espacio = "BP", qvalueCutoff = 0.05,
                       output_dir = "results/GO/pruned")
```

---

#### `graficar_go_balones()`

Visualizes a named list of GO enrichment results as a balloon/bubble chart where bubble size encodes fold enrichment and fill color encodes -log10(q-value).

| Parameter | Type | Description |
|---|---|---|
| `resuGO` | `list` | Named list of `enrichResult` objects (one per comparison) |

```r
p <- graficar_go_balones(go_pruned)
print(p)
```

---

### Save helpers

Three lightweight wrappers around `ggsave` that apply standardized dimensions for common plot types. All write into `output_dir`, which must be defined in the calling environment.

---

#### `save_pdf()`

Saves any ggplot (UMAP, FeaturePlot, etc.) as a PDF at 300 dpi with default dimensions of 10 x 8 inches.

| Parameter | Type | Description |
|---|---|---|
| `plot` | `ggplot` | Plot object to save |
| `file` | `character` | Filename appended to `output_dir` |
| `w` | `numeric` | Width in inches (default `10`) |
| `h` | `numeric` | Height in inches (default `8`) |

```r
save_pdf(umap_plot, "umap_clusters.pdf")
```

---

#### `save_vln()`

Saves a VlnPlot as a PDF via `save_pdf`, automatically scaling the height by the number of features plotted (14 x 6n inches).

| Parameter | Type | Description |
|---|---|---|
| `plot` | `ggplot` | VlnPlot object to save |
| `file` | `character` | Filename appended to `output_dir` |
| `n` | `integer` | Number of genes/features in the plot (default `1`) |

```r
save_vln(vln_plot, "vln_markers.pdf", n = 3)
```

---

#### `save_qc()`

Stacks a list of QC plots into a single column with `patchwork::wrap_plots` and saves the result as a PDF at 300 dpi with width 14 and height 6 per panel.

| Parameter | Type | Description |
|---|---|---|
| `plot_list` | `list` | List of ggplot objects to stack vertically |
| `file` | `character` | Filename appended to `output_dir` |

```r
save_qc(list(p_wt, p_mut), "qc_violin_grid.pdf")
```

---

## Quick start

```bash
# 1. Clone the repository
git clone https://github.com/your-org/ScRNASeq-Docker.git
cd ScRNASeq-Docker

# 2. Build the Docker image (R + all dependencies pre-installed)
docker build -t scrnaseq-docker .

# 3. Launch an interactive R session with the project directory mounted
docker run --rm -it \
  -v "$(pwd)":/workspace \
  -w /workspace \
  scrnaseq-docker R

# 4. Inside R: load all dependencies then source the function library
source("load_libraries.R")
source("ScRNA_Analysis_Functions.R")

# 5. Optionally source the custom Seurat extensions and run the full pipeline
source("custom_seurat.R")
source("scrnaseq_pipeline.R")
```

---

## Adapting to other organisms

Swap the organism-specific parameters below when applying this pipeline to species other than *Arabidopsis thaliana*. The `mt_pattern` and `cp_pattern` arguments are passed to `load_sample()` / `process_sample()`; the `orgdb` and `keytype` arguments are passed to `correr_enriquecimiento_go()`.

| Organism | `mt_pattern` | `cp_pattern` | `orgdb` | `keytype` |
|---|---|---|---|---|
| *Arabidopsis thaliana* | `"^ATMG"` | `"^ATCG"` | `org.At.tair.db` | `"TAIR"` |
| *Homo sapiens* | `"^MT-"` | `NULL` | `org.Hs.eg.db` | `"ENSEMBL"` |
| *Mus musculus* | `"^mt-"` | `NULL` | `org.Mm.eg.db` | `"ENSEMBL"` |

For rice (*Oryza sativa*) use `org.Os.eg.db` with `keytype = "GID"`, and set `mt_pattern` / `cp_pattern` to match your genome annotation's organelle gene naming convention. Also update the `universo` background gene vector in the GO enrichment step to reflect the full gene set of your target organism.


