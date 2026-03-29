# ScRNASeq Docker

Reproducible Docker environment for **Single-Cell RNA-seq analysis** with R 4.5 and Python 3.

---

## Table of Contents

- [Requirements](#requirements)
- [Quick Start](#quick-start)
- [Interactive Use](#interactive-use)
- [Included Packages](#included-packages)
- [Can I install other programs? (CellRanger, CellBender...)](#can-i-install-other-programs)
- [Build from Scratch](#build-from-scratch)
- [Repository Structure](#repository-structure)

---

## Requirements

- [Docker](https://docs.docker.com/get-docker/) installed
- ~5 GB of free disk space

---

## Quick Start

```bash
docker pull cliford2001/scrnaseq_docker:latest
```

---

## Interactive Use

**Interactive R session:**
```bash
docker run --rm -it -v $(pwd):/workspace cliford2001/scrnaseq_docker:latest
```

**Interactive Python session:**
```bash
docker run --rm -it -v $(pwd):/workspace cliford2001/scrnaseq_docker:latest python3
```

**Bash shell (explore the container):**
```bash
docker run --rm -it -v $(pwd):/workspace cliford2001/scrnaseq_docker:latest bash
```

**Run an R script:**
```bash
docker run --rm -v $(pwd):/workspace cliford2001/scrnaseq_docker:latest Rscript /workspace/my_script.R
```

**Run a Python script:**
```bash
docker run --rm -v $(pwd):/workspace cliford2001/scrnaseq_docker:latest python3 /workspace/my_script.py
```

> **Note:** The `-v $(pwd):/workspace` flag mounts your current folder inside the container at `/workspace`. Always run the command from the directory where your data is located.

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
| scater | 1.38.1 | QC and visualization (Bioconductor) |
| SingleCellExperiment | 1.32.0 | Bioconductor single-cell data structure |
| SummarizedExperiment | 1.40.0 | Bioconductor experiment container |
| zellkonverter | 1.20.1 | SCE ↔ AnnData conversion |
| harmony | 1.2.4 | Dataset integration / batch correction |
| DoubletFinder | 2.0.6 | Doublet detection |
| clustree | 0.5.1 | Clustering resolution visualization |
| hdf5r | 1.3.12 | HDF5 file I/O |
| Matrix | 1.7-5 | Sparse matrix support |
| ggplot2 | 4.0.2 | Data visualization |
| patchwork | 1.3.2 | Composing plots |
| cowplot | 1.2.0 | Publication-ready figures |
| gridExtra | 2.3 | Plot grids |
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

| Package | Description |
|---|---|
| scanpy | Single-cell analysis in Python |
| scFates | Trajectory analysis |
| palantir | Cell differentiation and trajectories |
| pandas | Data manipulation |
| numpy | Numerical computing |
| scipy | Scientific statistics |
| scikit-learn | Machine learning |
| matplotlib | Visualization |
| seaborn | Statistical visualization |
| cellbender | Ambient RNA removal |

---

## Can I install other programs?

**Yes.** Any tool available as a system package, pip package, or binary can be added to the Dockerfile.

### CellBender

CellBender (ambient RNA removal) is **already included** in this image. Use it directly:

```bash
cellbender remove-background \
    --input raw_feature_bc_matrix.h5 \
    --output output.h5
```

### CellRanger

CellRanger 10.0.0 is **already included** in this image. Use it directly:

```bash
docker run --rm -v $(pwd):/workspace scrnaseq_docker:latest \
    cellranger count \
    --id=sample \
    --transcriptome=/workspace/refdata \
    --fastqs=/workspace/fastqs \
    --sample=sample_name
```

### Other tools

| Tool | How to add |
|---|---|
| STARsolo | `apt-get install star` |
| Salmon / Alevin | `apt-get install salmon` |
| samtools | `apt-get install samtools` |
| FastQC | `apt-get install fastqc` |
| Trim Galore | `apt-get install trim-galore` |

---

## Build from Scratch

You need a [GitHub personal access token](https://github.com/settings/tokens) to install packages from GitHub.

```bash
git clone https://github.com/cliford2001/ScRNASeq-Docker.git
cd ScRNASeq-Docker
docker build --build-arg GITHUB_PAT=your_token -t scrnaseq_docker:latest .
```

> The build takes approximately **45–60 minutes** on the first run due to C++ compilation (BPCells/monocle3).

---

## Repository Structure

```
ScRNASeq-Docker/
├── Dockerfile          # Full image definition
├── docker-compose.yml  # docker compose configuration
└── workspace/          # Default mount point for your data
```

---

## Usage Examples

**R:**
```r
library(Seurat)
library(harmony)
library(monocle3)
library(DESeq2)
library(ggplot2)
library(patchwork)

# Load data (must be in /workspace)
seurat_obj <- readRDS("/workspace/my_object.rds")
```

**Python:**
```python
import scanpy as sc
import scFates as scf
import palantir
import pandas as pd
import numpy as np

# Load data
adata = sc.read_h5ad("/workspace/my_data.h5ad")
```

---

## Author

Developed by [@cliford2001](https://github.com/cliford2001)
