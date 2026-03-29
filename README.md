# ScRNASeq Docker

Reproducible Docker environment for **Single-Cell RNA-seq analysis** with R 4.5 and Python 3.

Includes: Seurat 5, monocle3, DESeq2, harmony, CellRanger 9.0.1, CellBender, scanpy, and more.

---

## Table of Contents

- [Requirements](#requirements)
- [Installing Docker](#installing-docker)
- [Quick Start](#quick-start)
- [Interactive Use](#interactive-use)
- [Included Packages](#included-packages)
- [CellRanger & CellBender](#cellranger--cellbender)
- [Build from Scratch](#build-from-scratch)
- [Repository Structure](#repository-structure)

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
| cellbender | Ambient RNA removal |
| pandas | Data manipulation |
| numpy | Numerical computing |
| scipy | Scientific statistics |
| scikit-learn | Machine learning |
| matplotlib | Visualization |
| seaborn | Statistical visualization |

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

## Repository Structure

```
ScRNASeq-Docker/
├── Dockerfile          # Full image definition
├── docker-compose.yml  # docker compose configuration
└── workspace/          # Default mount point for your data
```

---

## Author

Developed by [@cliford2001](https://github.com/cliford2001)
