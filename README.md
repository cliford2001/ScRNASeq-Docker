# ScRNASeq-Docker — Single-Cell RNA-Seq Analysis Pipeline

Helper library and Docker environment for single-cell RNA-seq of *Arabidopsis thaliana* (readily adaptable to other organisms). Covers the full trajectory from raw FASTQs to pseudotime.

---

## System requirements

### Minimum hardware

| Resource | Minimum | Recommended |
|---|---|---|
| RAM | 64 GB | 128 GB+ |
| CPU | 8 cores | 32+ cores |
| Storage | 100 GB free | 500 GB SSD |
| OS | Ubuntu 20.04 / macOS 13+ with Docker | Ubuntu 22.04+ |

> CellRanger's `--localcores` flag should match available CPUs. The pipeline was validated on a 24-core / 70 GB RAM server (R + Python steps) and a 128-core / 1 TB RAM server (CellRanger alignment with `--localcores=80`).

### Storage breakdown (4-sample run, Arabidopsis)

Measured from the actual pipeline run. Sizes may vary with sequencing depth and cell count.

| Component | Size |
|---|---|
| Docker image (`matigara/scrnaseq:latest`) | 10.4 GB |
| CellRanger 7.1.0 binary | 1.9 GB |
| Reference genome (TAIR10, `cellranger mkref` output) | 1.6 GB |
| Raw FASTQs (4 samples, R1 + R2) | 15 – 25 GB |
| CellRanger output (4 samples, `--no-bam`) | 2.5 GB |
| Seurat checkpoint objects (`.rds`) | 1.4 GB |
| Pipeline results (plots, tables, PDFs) | 2.2 GB |
| **Total** | **~35 – 45 GB** |

Key savings:
- `--no-bam` in `cellranger count` avoids BAM files (~20–30 GB per sample)
- Intermediate `.rds` checkpoints can be deleted once downstream steps are validated

---

## Repository structure

```
ScRNASeq-Docker/
├── Dockerfile                           # R 4.5.3 + Python 3.12.3 + CellRanger 7.1.0
├── docker-compose.yml                   # interactive session shortcut
├── README.md
├── README_capitulo3.md
└── workflow/
    ├── capitulo1_single_cell.R          # Chapter 1 — QC → clustering → annotation → export
    ├── capitulo2_pseudobulk_de.R        # Chapter 2 — DE → GO → hdWGCNA → TF network
    ├── capitulo3_pseudotime.ipynb       # Chapter 3 — trajectory → pseudotime (Python)
    ├── ScRNA_Analysis_Functions.R       # R helper function library (documented below)
    ├── ScRNA_Pseudotime_Functions.py    # Python helpers for Chapter 3
    ├── load_libraries.R                 # R package loader
    ├── load_libraries_python.py         # Python package loader
    └── custom_seurat.R                  # custom Seurat plot utilities
```

---

## Pipeline overview

```
Raw FASTQ files
       │
       ▼  Chapter 0 — CellRanger 7.1.0  (bash)
       │  mkref → count → filtered_feature_bc_matrix/
       │
       ▼  Chapter 1 — capitulo1_single_cell.R  (R / Seurat)
       │  QC → filtering → Harmony → clustering → annotation → .rds + .h5ad
       │
       ├─────────────────────────────────────────┐
       │                                         │
       ▼  Chapter 2 — capitulo2_pseudobulk_de.R  ▼  Chapter 3 — capitulo3_pseudotime.ipynb
       │  DESeq2 → volcano → GO enrichment        │  scFates trajectory → gene trends
       │  hdWGCNA co-expression network           │  pseudotime branches + milestones
       │  TF network (GENIE3 + WGCNA)             │
```

---

## Quick Start

### Option 1 — Docker Hub (recommended)

```bash
# Pull the pre-built image
docker pull matigara/scrnaseq:latest

# Launch an interactive bash session with your data mounted at /workspace
docker run -it --rm \
  -v /path/to/your/data:/workspace \
  matigara/scrnaseq:latest /bin/bash
```

Inside the container the pipeline scripts expect:

```r
PIPELINE_DIR <- "/workspace/ScRNASeq-Docker/workflow"
DATA_DIR     <- "/workspace/."
```

### Option 2 — docker compose

```bash
git clone https://github.com/cliford2001/ScRNASeq-Docker.git
cd ScRNASeq-Docker

docker compose run --rm r         # interactive R session
docker compose run --rm python    # interactive Python session
```

### Option 3 — Build locally

CellRanger requires manual download from [10x Genomics](https://www.10xgenomics.com/support/software/cell-ranger/downloads) (free registration). Place the tarball in the repo root before building:

```bash
# 1. Download cellranger-7.1.0.tar.gz from 10x Genomics and place it here
cp /downloads/cellranger-7.1.0.tar.gz .

# 2. Build the image
docker build -t scrnaseq:local .

# 3. Run
docker run -it --rm -v /path/to/your/data:/workspace scrnaseq:local /bin/bash
```

---

## Chapter 0 — CellRanger preprocessing (bash)

CellRanger 7.1.0 converts raw FASTQ files into the filtered feature-barcode matrices consumed by Chapter 1. Two steps are required: building a genome reference index and running per-sample alignment.

### Step 1 — Build the genome reference index

```bash
# ── Download genome + annotation (Arabidopsis thaliana TAIR10 — Ensembl Plants r59) ──
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-59/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-59/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.59.gtf.gz

gunzip Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
gunzip Arabidopsis_thaliana.TAIR10.59.gtf.gz

# ── Build the CellRanger reference (run once per genome build) ──
cellranger mkref \
  --genome=Arabidopsis_thaliana_TAIR10 \
  --fasta=Arabidopsis_thaliana.TAIR10.dna.toplevel.fa \
  --genes=Arabidopsis_thaliana.TAIR10.59.gtf \
  --nthreads=16
# Output: Arabidopsis_thaliana_TAIR10/   ← pass this path to --transcriptome below
```

> **Other organisms:** swap the FASTA/GTF for your target species. Update `mt_pattern`, `cp_pattern`, `orgdb`, and `keytype` in the R pipeline accordingly (see [Adapting to other organisms](#adapting-to-other-organisms)).

### Step 2 — Prepare FASTQs

CellRanger expects files named:

```
{SAMPLE}_S{N}_L001_R1_001.fastq.gz    # Read 1 — barcode + UMI  (28 bp)
{SAMPLE}_S{N}_L001_R2_001.fastq.gz    # Read 2 — cDNA insert
```

Rename if your files use a different convention:

```bash
# Example: public repository files → CellRanger convention
mv CRR775252_f1.fq.gz  scDS1a_S1_L001_R1_001.fastq.gz
mv CRR775252_r2.fq.gz  scDS1a_S1_L001_R2_001.fastq.gz
mv CRR775253_f1.fq.gz  scDS1b_S1_L001_R1_001.fastq.gz
mv CRR775253_r2.fq.gz  scDS1b_S1_L001_R2_001.fastq.gz
```

### Step 3 — Run alignment

```bash
#!/bin/bash
# run_cellranger.sh
set -euo pipefail

REF=/path/to/Arabidopsis_thaliana_TAIR10   # output of cellranger mkref
CR=cellranger                               # or /opt/cellranger-7.1.0/cellranger

# Skip samples whose output already exists
run_if_needed() {
  local id="$1"
  if [ -f "${id}/outs/metrics_summary.csv" ]; then
    echo "[$(date '+%F %T')] ${id}: already complete — skipping."
    return 0
  fi
  echo "[$(date '+%F %T')] Starting ${id} ..."
  "$CR" count \
    --id="${id}" \
    --fastqs=. \
    --sample="${id}" \
    --transcriptome="${REF}" \
    --localcores=80 \
    --no-bam
  echo "[$(date '+%F %T')] ${id}: done."
}

run_if_needed scDS1a
run_if_needed scDS1b
run_if_needed scDS2a
run_if_needed scDS2b
```

**Key flags:**

| Flag | Value | Description |
|---|---|---|
| `--id` | `scDS1a` | Name of the output directory |
| `--fastqs` | `.` | Directory containing FASTQ files |
| `--sample` | `scDS1a` | Prefix to match `{SAMPLE}_S*_L*_R*.fastq.gz` |
| `--transcriptome` | `/path/to/ref` | Directory produced by `cellranger mkref` |
| `--localcores` | `80` | CPU threads to use |
| `--no-bam` | — | Skip BAM output (~50 % disk saving) |

### Step 4 — Output structure

```
scDS1a/
└── outs/
    ├── filtered_feature_bc_matrix/      ← input to Chapter 1
    │   ├── barcodes.tsv.gz
    │   ├── features.tsv.gz
    │   └── matrix.mtx.gz
    ├── raw_feature_bc_matrix/
    ├── metrics_summary.csv              ← per-sample QC (cells detected, saturation)
    └── web_summary.html                 ← interactive QC report
```

The `filtered_feature_bc_matrix/` directory is loaded in Chapter 1 via:

```r
samples <- list(
  list(file = "cellranger_v2/scDS1a/outs/filtered_feature_bc_matrix",
       label = "scDS1a", condition = "condition_1"),
  list(file = "cellranger_v2/scDS1b/outs/filtered_feature_bc_matrix",
       label = "scDS1b", condition = "condition_1"),
  list(file = "cellranger_v2/scDS2a/outs/filtered_feature_bc_matrix",
       label = "scDS2a", condition = "condition_2"),
  list(file = "cellranger_v2/scDS2b/outs/filtered_feature_bc_matrix",
       label = "scDS2b", condition = "condition_2")
)
seurat_list_raw <- load_seurat_samples(samples, DATA_DIR,
                                       mt_pattern = "^ATMG",
                                       cp_pattern  = "^ATCG")
```

---

## Function reference

All functions live in `workflow/ScRNA_Analysis_Functions.R`. Source it at the top of any script:

```r
source(file.path(PIPELINE_DIR, "load_libraries.R"))
source(file.path(PIPELINE_DIR, "ScRNA_Analysis_Functions.R"))
```

The `save_pdf` / `save_vln` / `save_qc` helpers write to the global `output_dir` variable; reassign it before each section.

---

### 1 — Pipeline setup

| Function | Description |
|---|---|
| `create_pipeline_dirs(base_dir)` | Creates `01_qc/` … `09_pseudotime/` and `objects/` under `base_dir`; returns a named list of paths |

```r
list2env(create_pipeline_dirs(base_dir), envir = .GlobalEnv)
# creates: dir_01 … dir_09, dir_objects
```

---

### 2 — Data loading and QC

| Function | Description |
|---|---|
| `load_seurat_samples(samples, DATA_DIR, mt_pattern, cp_pattern)` | Batch-loads CellRanger matrices, computes organelle %, returns a named list of Seurat objects |
| `plot_qc_batch(seurat_list, colors, file)` | Saves a multi-panel QC violin grid (one panel per sample) to `output_dir/file` |
| `plot_qc_violin_grid(obj, label, color)` | Single-sample QC violin: nFeature, nCount, percent.mt, percent.cp |
| `filter_seurat_samples(seurat_list, min_features, max_mt)` | Applies thresholds + DoubletFinder to every sample in the list |

```r
# ── Standard workflow (CellRanger matrices) ──
seurat_list_raw <- load_seurat_samples(samples, DATA_DIR,
                                       mt_pattern = "^ATMG",
                                       cp_pattern  = "^ATCG")
plot_qc_batch(seurat_list_raw, colors, "qc_prefilter.pdf")

seurat_list <- filter_seurat_samples(seurat_list_raw,
                                     min_features = 200, max_mt = 5)
plot_qc_batch(seurat_list, colors, "qc_postfilter.pdf")
```

**For CellBender-filtered inputs** use the lower-level functions:

| Function | Description |
|---|---|
| `load_cellbender_filtered_h5(h5_path, project)` | Reads a CellBender `.h5` file into a Seurat object |
| `load_sample(sample_info, mt_pattern, cp_pattern)` | Load + annotate organelle %, no filtering |
| `filter_sample(obj, min_features, max_features, max_mt, max_cp, run_doubletfinder)` | Apply QC thresholds + optional DoubletFinder |
| `process_sample(sample_info, ...)` | `load_sample` + `filter_sample` in one call |

---

### 3 — Annotation

| Function | Description |
|---|---|
| `find_markers(seurat_obj, output_file, only_pos, min_pct, logfc_threshold, force)` | Runs `FindAllMarkers`; loads from TSV cache on subsequent calls unless `force = TRUE` |
| `annotate_by_markers(seurat_obj, markers, reference_file)` | Assigns cell types by crossing cluster markers with a tab-separated reference (`gene \| cell.types`); stores result in `$celltype` |
| `annotate_by_reference(seurat_obj, reference_obj, reference_col, dims)` | Label transfer via `FindTransferAnchors` + `TransferData`; result in `$celltype_reference` |
| `plot_marker_dotplot(seurat_obj, marker_table, annot_col, outfile, width, height)` | Dot plot of bibliography marker genes across clusters/cell types |

```r
markers <- find_markers(pbmc_harmony,
                        output_file = file.path(dir_03, "FindAllMarkers.tsv"))

pbmc_harmony <- annotate_by_markers(pbmc_harmony, markers,
                                    reference_file = "biblio_marks.txt")

pbmc_harmony <- annotate_by_reference(pbmc_harmony,
                                      reference_obj = ref_atlas,
                                      reference_col = "annotation")

plot_marker_dotplot(pbmc_harmony, marker_table,
                    annot_col = "celltype",
                    outfile   = file.path(dir_03, "dotplot_biblio.pdf"))
```

---

### 4 — Curation and subclustering

| Function | Description |
|---|---|
| `subcluster_cell_type(obj, cell_type, annot_col, resolution, dims)` | Subsets to one cell type and re-runs PCA → UMAP → clustering |
| `plot_subcluster_umap(subcluster_obj, cell_type, output_dir)` | Saves and returns a UMAP colored by sub-cluster |
| `save_subcluster_composite(subcluster_list, marker_table, output_dir)` | Composite PDF with [UMAP \| marker dotplot] per cell type |
| `apply_subcluster_reassignment(obj, subcluster_list, reassign, source_col, dest_col)` | Applies a sub-cluster → cell-type reassignment map to the global object |

```r
mesophyll_sub <- subcluster_cell_type(pbmc_harmony, "Mesophyll",
                                      annot_col = "celltype_grouped")
p_meso <- plot_subcluster_umap(mesophyll_sub, "Mesophyll", dir_05)

pbmc_harmony <- apply_subcluster_reassignment(
  obj             = pbmc_harmony,
  subcluster_list = list(mesophyll_sub = mesophyll_sub),
  reassign        = list(mesophyll_sub = c("0" = "Mesophyll", "others" = "Mesophyll")),
  source_col      = "celltype_grouped",
  dest_col        = "celltype_curated"
)
```

---

### 5 — Export

| Function | Description |
|---|---|
| `export_to_scanpy(seurat_obj, outfile)` | Converts Seurat → SingleCellExperiment → `.h5ad` (zellkonverter); embeds PCA, UMAP, and Harmony reductions |

```r
export_to_scanpy(pbmc_harmony,
                 file.path(dir_objects, "pbmc_harmony_curated.h5ad"))
```

---

### 6 — Pseudobulk, DESeq2, and visualization

| Function | Description |
|---|---|
| `create_cell_type_subsets(seurat_obj, annot_col)` | Splits object into per-cell-type named list |
| `assign_pseudoreplicates_batch(cell_type_subsets, n_reps, conditions, seed)` | Assigns pseudo-replicate labels within each condition per cell type |
| `run_pseudobulk_deseq2_analysis(cell_type_subsets_replicates, comparisons, output_dir)` | Aggregates counts, runs DESeq2, writes per-comparison CSVs |
| `render_volcano_plots(results_dir, padj_cut, lfc_cut, output_dir)` | Batch volcano plots for every DESeq2 result file |
| `build_differential_tables(results_dir, padj_cut, lfc_cut, output_dir)` | Filtered significant-gene tables with up/down classification |
| `build_logfc_heatmap(logfc_table, output_dir, ...)` | Hierarchically clustered log2FC heatmap across cell types |
| `plot_volcano(file, padj_cut, lfc_cut)` | Single volcano plot from a DESeq2 CSV |
| `plot_heatmap(matriz, min_genes, deepSplit_val, breaks)` | Dynamic-tree-cut heatmap |
| `plot_replicate_correlation(pseudobulk_mat, main)` | Pearson correlation heatmap across pseudo-replicates |

```r
cell_type_subsets <- create_cell_type_subsets(pbmc_harmony,
                                               annot_col = "celltype_grouped")

cell_type_subsets_replicates <- assign_pseudoreplicates_batch(
  cell_type_subsets, n_reps = 3,
  conditions = c("condition_1", "condition_2"), seed = 1807)

run_pseudobulk_deseq2_analysis(
  cell_type_subsets_replicates,
  comparisons = list(list(conds = c("condition_1", "condition_2"),
                          tag   = "cond1_vs_cond2")),
  output_dir = dir_06)

render_volcano_plots(dir_06, padj_cut = 0.05, lfc_cut = 1,
                     output_dir = dir_06)
```

---

### 7 — GO enrichment

| Function | Description |
|---|---|
| `run_simple_go_enrichment(diff_table, universe, orgdb, keytype, output_dir)` | Per-cell-type GO enrichment from a differential gene table |
| `run_go_enrichment_suite(diff_table, universe, orgdb, keytype, namespaces, output_dir)` | Full GO suite across BP/MF/CC with optional pruning |
| `correr_enriquecimiento_go(tabla, universo, espacio, orgdb, keytype, ...)` | Low-level: enrichGO on a binary gene matrix |
| `podar_go(resuGO, nivel, ...)` | Filter enrichResult list by GO level |
| `graficar_go_balones(resuGO)` | Bubble chart: fold-enrichment × −log10(q-value) |

```r
run_simple_go_enrichment(
  diff_table = de_results,
  universe   = rownames(pbmc_harmony),
  orgdb      = org.At.tair.db,
  keytype    = "TAIR",
  output_dir = dir_07)
```

---

### 8 — hdWGCNA co-expression network

| Function | Description |
|---|---|
| `run_hdwgcna(seurat_obj, output_dir, soft_power, ...)` | Full hdWGCNA pipeline: metacell construction → module detection → ME correlation |
| `plot_hdwgcna_network(hdwgcna_dir, output_dir, tom_threshold, max_modules, ...)` | ggraph network plots per module, filtered by TOM weight and DE genes |
| `filter_hdwgcna_by_de(hdwgcna_dir, de_dir, output_dir, tom_threshold, ...)` | Prunes network to DE genes only and re-plots |

```r
run_hdwgcna(pbmc_harmony,
            output_dir  = dir_08,
            soft_power  = NULL)   # auto-detect soft power (R² ≥ 0.8)

plot_hdwgcna_network(hdwgcna_dir = dir_08,
                     output_dir  = file.path(dir_08, "network_wgcna"),
                     tom_threshold = 0.2)

filter_hdwgcna_by_de(hdwgcna_dir = dir_08,
                     de_dir      = dir_06,
                     output_dir  = file.path(dir_08, "network_wgcna_DE"),
                     tom_threshold = 0.2)
```

---

### 9 — Save helpers

All three write into `output_dir` (must be defined in the calling scope before use).

| Function | Signature | Description |
|---|---|---|
| `save_pdf(plot, file, w, h)` | `w=10, h=8` | Standard plot — UMAP, FeaturePlot, etc. |
| `save_vln(plot, file, n)` | `n=1` | Violin plot; `n` = number of genes plotted |
| `save_qc(plot_list, file)` | — | Stacks a list of plots vertically |

```r
output_dir <- dir_03   # reassign before each section

save_pdf(DimPlot(pbmc_harmony, group.by = "celltype", label = TRUE),
         "umap_annotation_biblio.pdf")

save_vln(VlnPlot(pbmc_harmony, features = c("AT5G26000", "AT5G54250")),
         "vln_guard_genes.pdf", n = 2)
```

---

## Adapting to other organisms

Change the three organism-specific arguments across all chapters:

| Organism | `mt_pattern` | `cp_pattern` | `orgdb` | `keytype` |
|---|---|---|---|---|
| *Arabidopsis thaliana* | `"^ATMG"` | `"^ATCG"` | `org.At.tair.db` | `"TAIR"` |
| *Homo sapiens* | `"^MT-"` | `NULL` | `org.Hs.eg.db` | `"ENSEMBL"` |
| *Mus musculus* | `"^mt-"` | `NULL` | `org.Mm.eg.db` | `"ENSEMBL"` |
| *Oryza sativa* | `"^LoChO"` | `NULL` | `org.Os.eg.db` | `"GID"` |

Also update the `universe` background gene vector in `run_simple_go_enrichment` to the full gene set of your organism.

---

## Software versions

All analyses run inside `matigara/scrnaseq:latest` (Ubuntu 24.04.4 LTS, R 4.5.3, Python 3.12.3, CellRanger 7.1.0). CellRanger is bundled in the image — no separate installation required.

### R packages

| Package | Version | Package | Version |
|---|---|---|---|
| Seurat | 5.5.0 | SeuratObject | 5.4.0 |
| SeuratDisk | 0.0.0.9021 | harmony | 2.0.3 |
| hdWGCNA | 0.4.11 | DESeq2 | 1.50.2 |
| clusterProfiler | 4.18.4 | org.At.tair.db | 3.22.0 |
| DoubletFinder | 2.0.6 | clustree | 0.5.1 |
| scater | 1.38.1 | SingleCellExperiment | 1.32.0 |
| zellkonverter | 1.20.1 | ggplot2 | 4.0.3 |
| ggrepel | 0.9.8 | ggraph | 2.2.2 |
| patchwork | 1.3.2 | cowplot | 1.2.0 |
| dplyr | 1.2.1 | tidyverse | 2.0.0 |
| igraph | 2.3.2 | tidygraph | 1.3.1 |
| WGCNA | 1.74 | dynamicTreeCut | 1.63.1 |
| enrichR | 3.4 | UCell | 2.14.0 |
| Matrix | 1.7.5 | hdf5r | 1.3.12 |
| reticulate | 1.46.0 | BiocGenerics | 0.56.0 |

### Python packages

| Package | Version | Package | Version |
|---|---|---|---|
| scanpy | 1.12.1 | anndata | 0.12.16 |
| scFates | 1.2.4 | palantir | 1.4.4 |
| numpy | 2.4.6 | pandas | 2.3.3 |
| scipy | 1.17.1 | matplotlib | 3.10.9 |
| seaborn | 0.13.2 | scikit-learn | 1.9.0 |
| umap-learn | 0.5.12 | igraph | 1.0.0 |
