# ScRNASeq Docker

Reproducible Docker environment for **Single-Cell RNA-seq analysis** with R 4.5 and Python 3.

This pipeline is **organism-agnostic**: it ships with defaults tuned for *Arabidopsis thaliana*, but running it on human, mouse, or any other species requires changing only a few parameters (see [Adapting to Other Organisms](#adapting-to-other-organisms)). The full workflow covers raw-count loading (CellBender output), QC, doublet detection, batch correction, clustering, cell-type annotation, pseudobulk differential expression, volcano plots, hierarchical heatmaps, and GO enrichment — all driven from a single main script (`paperrr.R`) and a reusable functions library (`ScRNA_Analysis_Functions.R`).

Includes: Seurat 5, monocle3, DESeq2, harmony, CellRanger 9.0.1, CellBender, scanpy, and more.

---

## Table of Contents

- [Requirements](#requirements)
- [Installing Docker](#installing-docker)
- [Quick Start](#quick-start)
- [Interactive Use](#interactive-use)
- [Pipeline Overview](#pipeline-overview)
- [Adapting to Other Organisms](#adapting-to-other-organisms)
- [Included Packages](#included-packages)
- [CellRanger & CellBender](#cellranger--cellbender)
- [Build from Scratch](#build-from-scratch)
- [Repository Structure](#repository-structure)
- [Troubleshooting](#troubleshooting)
- [Author](#author)

---

## Requirements

- Docker installed (see below)
- ~10 GB of free disk space

---

## Installing Docker

### Linux (Ubuntu/Debian)

```bash
sudo apt-get update
sudo apt-get install -y docker.io
sudo systemctl enable --now docker
sudo usermod -aG docker $USER   # allows running Docker without sudo
```
> **Re-login** after running `usermod` for the change to take effect.
> Verify: `docker --version`

### Mac

1. Download [Docker Desktop for Mac](https://www.docker.com/products/docker-desktop/)
2. Open the `.dmg` and drag Docker to Applications
3. Launch Docker Desktop from Applications
4. Verify: open Terminal and run `docker --version`

### Windows (PowerShell)

1. Download [Docker Desktop for Windows](https://www.docker.com/products/docker-desktop/)
2. Run the installer (requires WSL 2 — the installer will guide you)
3. Launch Docker Desktop
4. Verify: open PowerShell and run:
```powershell
docker --version
```

---

## Quick Start

### Linux / Mac

```bash
docker pull mvergara19/scrnaseq_docker:latest
```

### Windows (PowerShell)

```powershell
docker pull mvergara19/scrnaseq_docker:latest
```

---

## Interactive Use

### Linux / Mac

**Interactive R session:**
```bash
docker run --rm -it -v $(pwd):/workspace mvergara19/scrnaseq_docker:latest
```

**Interactive Python session:**
```bash
docker run --rm -it -v $(pwd):/workspace mvergara19/scrnaseq_docker:latest python3
```

**Bash shell:**
```bash
docker run --rm -it -v $(pwd):/workspace mvergara19/scrnaseq_docker:latest bash
```

**Run an R script:**
```bash
docker run --rm -v $(pwd):/workspace mvergara19/scrnaseq_docker:latest Rscript /workspace/my_script.R
```

**Run a Python script:**
```bash
docker run --rm -v $(pwd):/workspace mvergara19/scrnaseq_docker:latest python3 /workspace/my_script.py
```

### Windows (PowerShell)

**Interactive R session:**
```powershell
docker run --rm -it -v ${PWD}:/workspace mvergara19/scrnaseq_docker:latest
```

**Interactive Python session:**
```powershell
docker run --rm -it -v ${PWD}:/workspace mvergara19/scrnaseq_docker:latest python3
```

**Bash shell:**
```powershell
docker run --rm -it -v ${PWD}:/workspace mvergara19/scrnaseq_docker:latest bash
```

**Run an R script:**
```powershell
docker run --rm -v ${PWD}:/workspace mvergara19/scrnaseq_docker:latest Rscript /workspace/my_script.R
```

**Run a Python script:**
```powershell
docker run --rm -v ${PWD}:/workspace mvergara19/scrnaseq_docker:latest python3 /workspace/my_script.py
```

> **Note:** The `-v` flag mounts your current folder into the container at `/workspace`.
> Linux/Mac use `$(pwd)`, Windows PowerShell uses `${PWD}`.
> Always run the command from the directory where your data is located.

---

## Pipeline Overview

The R pipeline (`paperrr.R`) runs the following steps in order:

1. **QC pre-filter** — violin plots of nFeature, nCount, percent.mt, and percent.cp per sample
2. **QC filtering** — remove low-quality cells based on minimum feature count and maximum mitochondrial %
3. **Merge & preprocessing** — normalize (log), find variable features, scale data, PCA, initial UMAP
4. **Harmony batch correction** — integrate all samples by `orig.ident` to remove batch effects
5. **Elbow plot** — k-means WSS curve to guide the choice of significant PCA dimensions
6. **Clustree resolution sweep** — run `FindClusters` across a range of resolutions and visualize stability
7. **Final clustering** — apply the chosen resolution (default 0.3) with the Leiden algorithm
8. **Cell count plots** — barplots of cells per cluster broken down by sample and condition
9. **Global pseudobulk** — aggregate counts across all cells per sample; Pearson replicate-correlation heatmap
10. **Marker genes** — `FindAllMarkers` with bibliography-based cell-type annotation
11. **Reference annotation** — transfer cell-type labels from an external reference Seurat object
12. **Annotated clustree** — resolution sweep overlaid with the transferred cell-type labels
13. **Gene expression visualisation** — VlnPlot and FeaturePlot for individual genes and gene sets, both globally and within a selected cell type
14. **Cell type grouping** — merge fine-grained annotations into broader categories via a configurable mapping table
15. **Interactive curation** — manual subclustering and reassignment of heterogeneous cell types (run interactively, not sourced)
16. **Export to h5ad** — convert the curated Seurat object to AnnData format for downstream Python/scanpy analysis
17. **Dotplot** — marker-gene dot plot across all curated cell types
18. **Pseudoreplicates & pseudobulk per cell type** — assign pseudo-replicates and aggregate counts; one CSV per cell type
19. **DESeq2** — all pairwise condition comparisons, one results subdirectory per comparison
20. **Volcano plots** — publication-ready volcano plots, two panels per page, one PDF per comparison
21. **Differential table & heatmap** — log2FC matrix across all cell types; hierarchical clustering heatmap
22. **GO enrichment** — clusterProfiler `enrichGO`, dot plots for raw and simplified terms, pruned GO trees

All major plots are saved as **PDF files** under `results/` or the configured `output_dir`.

---

## Adapting to Other Organisms

Three parameters control the organism-specific parts of the pipeline.

### Mitochondrial / chloroplast gene patterns (`load_sample`)

```r
# Arabidopsis thaliana (default)
mt_pattern = "^ATMG"
cp_pattern  = "^ATCG"

# Human
mt_pattern = "^MT-"
cp_pattern  = NULL

# Mouse
mt_pattern = "^mt-"
cp_pattern  = NULL
```

### GO enrichment database (`paperrr.R` — GO ENRICHMENT section)

```r
# Arabidopsis thaliana (default)
orgdb   <- org.At.tair.db
keytype <- "TAIR"

# Human
orgdb   <- org.Hs.eg.db
keytype <- "ENSEMBL"   # or "SYMBOL"

# Mouse
orgdb   <- org.Mm.eg.db
keytype <- "ENSEMBL"   # or "SYMBOL"
```

Also update the `universo` background gene list to match your organism.

### DESeq2 comparisons

The `comparaciones` list in `paperrr.R` is fully configurable. Add or remove entries to match your condition names and the contrasts you need:

```r
comparaciones <- list(
  list(conds = c("Control", "Treatment"), tag = "ctrl_vs_treat"),
  list(conds = c("Control", "HighDose"),  tag = "ctrl_vs_high")
)
```

Each `tag` becomes a subdirectory under `results/deseq2/`, `results/diff/`, and `results/Enrichment/`.

---

## Included Packages

### R 4.5

| Package | Version | Description |
|---|---|---|
| Seurat | 5.4.0 | Single-cell analysis framework |
| SeuratDisk | 0.0.0.9021 | Read/write h5seurat and h5ad formats |
| SeuratWrappers | 0.4.0 | Integrations with external methods |
| monocle3 | 1.4.26 | Trajectory analysis |
| DESeq2 | 1.50.2 | Differential expression |
| clusterProfiler | — | GO and pathway enrichment |
| org.At.tair.db | — | Arabidopsis thaliana annotation database (swap for org.Hs.eg.db, org.Mm.eg.db, etc.) |
| scater | 1.38.1 | QC and visualization (Bioconductor) |
| SingleCellExperiment | 1.32.0 | Bioconductor single-cell data structure |
| SummarizedExperiment | 1.40.0 | Bioconductor experiment container |
| zellkonverter | 1.20.1 | SCE <-> AnnData conversion |
| harmony | 1.2.4 | Dataset integration / batch correction |
| DoubletFinder | 2.0.6 | Doublet detection |
| clustree | 0.5.1 | Clustering resolution visualization |
| dynamicTreeCut | — | Dynamic dendrogram-based cluster cutting |
| hdf5r | 1.3.12 | HDF5 file I/O |
| Matrix | 1.7-5 | Sparse matrix support |
| ggplot2 | 4.0.2 | Data visualization |
| ggrepel | — | Non-overlapping text labels |
| ggpubr | — | Publication-ready ggplot2 figures |
| pheatmap | — | Pretty heatmaps |
| RColorBrewer | — | Color palettes |
| patchwork | 1.3.2 | Composing plots |
| cowplot | 1.2.0 | Publication-ready figures |
| gridExtra | 2.3 | Plot grids |
| reshape2 | — | Data reshaping (melt/cast) |
| dplyr | 1.2.0 | Data manipulation |
| tibble | 3.3.1 | Modern data frames |
| tidyverse | 2.0.0 | Data analysis ecosystem |
| knitr | 1.51 | Dynamic reports |
| kableExtra | 1.4.0 | Enhanced tables in reports |
| VennDiagram | 1.8.2 | Venn diagrams |
| ggvenn | 0.1.19 | Venn diagrams with ggplot2 |
| eulerr | 7.0.4 | Euler diagrams |
| UpSetR | 1.4.0 | UpSet plots |

### Python 3

| Package | Version | Description |
|---|---|---|
| scanpy | 1.12 | Single-cell analysis in Python |
| scFates | 1.2.3 | Trajectory analysis |
| palantir | 1.4.4 | Cell differentiation and trajectories |
| cellbender | 0.3.0 | Ambient RNA removal |
| pandas | 2.3.3 | Data manipulation |
| numpy | 2.4.3 | Numerical computing |
| scipy | 1.17.1 | Scientific statistics |
| scikit-learn | 1.8.0 | Machine learning |
| matplotlib | 3.10.8 | Visualization |
| seaborn | 0.13.2 | Statistical visualization |

### Command-line tools

| Tool | Version | Description |
|---|---|---|
| CellRanger | 9.0.1 | 10x Genomics read alignment and quantification |

---

## CellRanger & CellBender

### CellRanger

CellRanger 9.0.1 is already included. Example usage:

**Linux / Mac:**
```bash
docker run --rm -v $(pwd):/workspace mvergara19/scrnaseq_docker:latest \
    cellranger count --id=sample_id --fastqs=/workspace/fastqs --sample=sample_name --transcriptome=/workspace/refdata
```

**Windows (PowerShell):**
```powershell
docker run --rm -v ${PWD}:/workspace mvergara19/scrnaseq_docker:latest `
    cellranger count --id=sample_id --fastqs=/workspace/fastqs --sample=sample_name --transcriptome=/workspace/refdata
```

### CellBender

CellBender is already included. Example usage:

**Linux / Mac:**
```bash
docker run --rm -v $(pwd):/workspace mvergara19/scrnaseq_docker:latest \
    cellbender remove-background --input /workspace/raw_feature_bc_matrix.h5 --output /workspace/output.h5
```

**Windows (PowerShell):**
```powershell
docker run --rm -v ${PWD}:/workspace mvergara19/scrnaseq_docker:latest `
    cellbender remove-background --input /workspace/raw_feature_bc_matrix.h5 --output /workspace/output.h5
```

---

## Build from Scratch

You need a [GitHub personal access token](https://github.com/settings/tokens) to install packages from GitHub.
You also need the `cellranger-9.0.1.tar.gz` file downloaded from [10x Genomics](https://www.10xgenomics.com/support/software/cell-ranger/downloads) placed in the same folder as the Dockerfile.

**Linux / Mac:**
```bash
git clone https://github.com/cliford2001/ScRNASeq-Docker.git
cd ScRNASeq-Docker
# Place cellranger-9.0.1.tar.gz here
docker build --build-arg GITHUB_PAT=your_token -t scrnaseq_docker:latest .
```

**Windows (PowerShell):**
```powershell
git clone https://github.com/cliford2001/ScRNASeq-Docker.git
cd ScRNASeq-Docker
# Place cellranger-9.0.1.tar.gz here
docker build --build-arg GITHUB_PAT=your_token -t scrnaseq_docker:latest .
```

> The build takes approximately **45–60 minutes** on the first run due to C++ compilation (BPCells/monocle3).

---

## Repository Structure

```
ScRNASeq-Docker/
├── Dockerfile                       # Full image definition
├── docker-compose.yml               # docker compose configuration
├── load_libraries.R                 # All library() calls (sourced by paperrr.R)
├── paperrr.R                        # Main analysis pipeline (source this script)
├── ScRNA_Analysis_Functions.R       # Reusable helper functions
├── custom_seurat.R                  # Custom Seurat utilities
└── workspace/                       # Default mount point for your data
```

### Pipeline output layout

```
results/
├── objs/                            # Exported h5ad objects (Seurat -> AnnData)
├── pseudobulk_replicas/             # Pseudobulk count matrices per cell type (.csv)
├── deseq2/
│   ├── 05_5/                        # DESeq2 results: 0.5N vs 5N
│   ├── 0_5/                         # DESeq2 results: 0N vs 5N
│   └── 0_05/                        # DESeq2 results: 0N vs 0.5N
├── volcano_plots/                   # Volcano plot PDFs, one per comparison
├── diff/
│   ├── 05_5/                        # Differential gene tables and heatmaps
│   ├── 0_5/
│   └── 0_05/
└── Enrichment/
    ├── 05_5/                        # GO enrichment dot plots
    ├── 0_5/
    └── 0_05/
```

Major QC and annotation PDFs (e.g. `qc_prefilter.pdf`, `umap_postharmony.pdf`, `clustree.pdf`, `dotplot_marcadores.pdf`) are written to the `output_dir` configured at the top of `paperrr.R`.

---

## Troubleshooting

### `permission denied` when running `docker pull` or `docker run`

**Error:**
```
permission denied while trying to connect to the Docker daemon socket at
unix:///var/run/docker.sock
```

**Cause:** Your user is not in the `docker` group. By default, only `root` can
communicate with the Docker daemon.

**Fix:** Run this once (requires admin/sudo on that machine):

```bash
sudo usermod -aG docker $USER
```

Then **log out and log back in** for the change to take effect. After that,
Docker works without `sudo` for that user permanently.

If you cannot log out (e.g., remote server), use this as a temporary workaround
for the current session only:

```bash
newgrp docker
```

> **Note:** This is an OS-level permission issue, not a problem with the Docker
> image. The image is public and requires no authentication to pull.

### Seurat 5: `slot` argument error in `GetAssayData` / `AggregateExpression`

**Error (example):**
```
Error in GetAssayData(..., slot = "counts") :
  'slot' is deprecated; use 'layer' instead
```

**Cause:** Seurat 5 replaced the `slot` argument with `layer` across its API.

**Fix:** This has already been corrected in `ScRNA_Analysis_Functions.R`.
Both `exportar_para_scanpy` and `hacer_pseudobulk` now use `layer = "counts"` /
`layer = "data"` instead of the old `slot` argument. If you encounter this error
in custom code, apply the same substitution.

---

## Author

Developed by [@cliford2001](https://github.com/cliford2001)
