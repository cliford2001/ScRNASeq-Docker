# ScRNASeq Docker

Entorno Docker reproducible para análisis de **Single-Cell RNA-seq** con R 4.5 y Python 3.

---

## Contenido

- [Requisitos](#requisitos)
- [Instalación rápida](#instalación-rápida)
- [Uso interactivo](#uso-interactivo)
- [Paquetes incluidos](#paquetes-incluidos)
- [Construir la imagen desde cero](#construir-la-imagen-desde-cero)
- [Estructura del repositorio](#estructura-del-repositorio)

---

## Requisitos

- [Docker](https://docs.docker.com/get-docker/) instalado
- ~5 GB de espacio en disco

---

## Instalación rápida

```bash
docker pull cliford2001/scrnaseq_docker:latest
```

---

## Uso interactivo

**Sesión R:**
```bash
docker run --rm -it -v $(pwd):/workspace cliford2001/scrnaseq_docker:latest
```

**Sesión Python:**
```bash
docker run --rm -it -v $(pwd):/workspace cliford2001/scrnaseq_docker:latest python3
```

**Bash (explorar el contenedor):**
```bash
docker run --rm -it -v $(pwd):/workspace cliford2001/scrnaseq_docker:latest bash
```

**Correr un script R:**
```bash
docker run --rm -v $(pwd):/workspace cliford2001/scrnaseq_docker:latest Rscript /workspace/mi_script.R
```

**Correr un script Python:**
```bash
docker run --rm -v $(pwd):/workspace cliford2001/scrnaseq_docker:latest python3 /workspace/mi_script.py
```

> **Nota:** El flag `-v $(pwd):/workspace` monta tu carpeta actual dentro del contenedor en `/workspace`. Siempre ejecuta el comando desde la carpeta donde están tus datos.

---

## Paquetes incluidos

### R 4.5

| Paquete | Versión | Descripción |
|---|---|---|
| Seurat | 5.4.0 | Análisis single-cell |
| SeuratDisk | 0.0.0.9021 | Lectura/escritura de formatos h5seurat/h5ad |
| SeuratWrappers | 0.4.0 | Integraciones con otros métodos |
| monocle3 | 1.4.26 | Análisis de trayectorias |
| DESeq2 | 1.50.2 | Expresión diferencial |
| scater | 1.38.1 | Control de calidad y visualización |
| SingleCellExperiment | 1.32.0 | Estructura de datos Bioconductor |
| SummarizedExperiment | 1.40.0 | Contenedor de experimentos |
| zellkonverter | 1.20.1 | Conversión SCE ↔ AnnData |
| harmony | 1.2.4 | Integración de datasets |
| DoubletFinder | 2.0.6 | Detección de dobletes |
| clustree | 0.5.1 | Visualización de resoluciones de clustering |
| hdf5r | 1.3.12 | Lectura/escritura HDF5 |
| Matrix | 1.7-5 | Matrices dispersas |
| ggplot2 | 4.0.2 | Visualización |
| patchwork | 1.3.2 | Composición de gráficos |
| cowplot | 1.2.0 | Figuras científicas |
| gridExtra | 2.3 | Grids de gráficos |
| dplyr | 1.2.0 | Manipulación de datos |
| tibble | 3.3.1 | Tablas modernas |
| tidyverse | 2.0.0 | Ecosistema de análisis de datos |
| knitr | 1.51 | Reportes dinámicos |
| kableExtra | 1.4.0 | Tablas en reportes |
| VennDiagram | 1.8.2 | Diagramas de Venn |
| ggvenn | 0.1.19 | Diagramas de Venn con ggplot2 |
| eulerr | 7.0.4 | Diagramas de Euler |
| UpSetR | 1.4.0 | Gráficos UpSet |
| harmony | 1.2.4 | Integración batch |

### Python 3

| Paquete | Descripción |
|---|---|
| scanpy | Análisis single-cell en Python |
| scFates | Análisis de trayectorias |
| palantir | Diferenciación celular y trayectorias |
| pandas | Manipulación de datos |
| numpy | Cálculo numérico |
| scipy | Estadística científica |
| scikit-learn | Machine learning |
| matplotlib | Visualización |
| seaborn | Visualización estadística |

---

## Construir la imagen desde cero

Necesitas un [token de GitHub](https://github.com/settings/tokens) para instalar los paquetes privados/GitHub.

```bash
git clone https://github.com/cliford2001/ScRNASeq-Docker.git
cd ScRNASeq-Docker
docker build --build-arg GITHUB_PAT=tu_token -t scrnaseq_docker:latest .
```

> El build tarda aproximadamente **45-60 minutos** la primera vez por la compilación de paquetes C++ (BPCells/monocle3).

---

## Estructura del repositorio

```
ScRNASeq-Docker/
├── Dockerfile          # Definición completa de la imagen
├── docker-compose.yml  # Configuración para docker compose
└── workspace/          # Carpeta montada por defecto (tus datos van aquí)
```

---

## Ejemplo de uso en R

```r
library(Seurat)
library(harmony)
library(monocle3)
library(DESeq2)
library(ggplot2)
library(patchwork)

# Cargar datos (deben estar en /workspace)
seurat_obj <- readRDS("/workspace/mi_objeto.rds")
```

## Ejemplo de uso en Python

```python
import scanpy as sc
import scFates as scf
import palantir
import pandas as pd
import numpy as np

# Cargar datos
adata = sc.read_h5ad("/workspace/mi_datos.h5ad")
```

---

## Autor

Desarrollado por [@cliford2001](https://github.com/cliford2001)
