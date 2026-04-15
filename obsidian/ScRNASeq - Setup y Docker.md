---
tags: [scrnaseq, docker, setup]
type: note
relacionado: [[ScRNASeq-Docker MOC]]
---

# ScRNASeq — Setup y Docker

## Imagen Docker

**Base**: `rocker/r-ver:4.5`  
**Nombre imagen**: `scrnaseq_docker:latest`  
**Nombre contenedor**: `r45`  
**Directorio de trabajo interno**: `/workspace`

---

## Cómo correr el contenedor

```bash
# Clonar el repo
git clone https://github.com/cliford2001/ScRNASeq-Docker.git
cd ScRNASeq-Docker

# Construir y levantar el contenedor interactivo
docker-compose up -d
docker exec -it r45 R
```

O directamente:
```bash
docker build -t scrnaseq_docker .
docker run -it --rm \
  -v $(pwd)/workspace:/workspace \
  scrnaseq_docker R
```

---

## Estructura del docker-compose.yml

```yaml
services:
  r:
    build: .
    image: scrnaseq_docker:latest
    container_name: r45
    volumes:
      - ./workspace:/workspace          # Proyecto local montado
      - r-packages:/usr/local/lib/R/library  # Cache de paquetes R persistente
    stdin_open: true
    tty: true
    command: R

volumes:
  r45-packages:    # Persiste paquetes entre reinicios del contenedor
```

> **Clave**: Los paquetes R instalados dentro del contenedor persisten en un volumen nombrado (`r45-packages`). No se reinstalan en cada `docker-compose up`.

---

## Dependencias del sistema (apt-get)

| Categoría | Paquetes |
|---|---|
| Compiladores | git, cmake, make, wget, curl |
| Criptografía/HTTP | libcurl4-openssl-dev, libssl-dev, libxml2-dev |
| Gráficos | libfontconfig1-dev, libharfbuzz-dev, libfribidi-dev, libfreetype6-dev, libpng-dev, libtiff5-dev, libjpeg-dev |
| Geoespacial | libgdal-dev, libgeos-dev, libproj-dev, libsqlite3-dev, libudunits2-dev |
| Científico | libhdf5-dev (HDF5), libv8-dev, libglpk-dev, libfftw3-dev, libgsl-dev |
| Git/SSH | libgit2-dev, libssh2-1-dev |
| Python | python3, python3-pip, python3-venv |

---

## Entorno Python (venv en /opt/venv)

```
scanpy          # Análisis single-cell en Python
scFates         # Inferencia de trayectorias
palantir        # Inferencia de destino celular
cellbender      # Remoción de RNA ambiente
pandas, numpy, scipy, scikit-learn
matplotlib, seaborn
```

---

## Paquetes R instalados

### CRAN (41 paquetes clave)

| Categoría | Paquetes |
|---|---|
| Single-cell | Seurat, SeuratObject, harmony |
| Visualización | ggplot2, patchwork, cowplot, gridExtra, ggrepel, ggbeeswarm, ggridges, ggforce, pheatmap |
| Datos | dplyr, tidyverse, tibble, data.table, reshape2 |
| Clustering | clustree, ggraph, igraph, dynamicTreeCut |
| I/O | hdf5r, Matrix |
| Paralelización | future, future.apply, parallelly |
| Interactividad | plotly, shiny, miniUI |
| Puente Python | reticulate |
| Reportes | rmarkdown, knitr |

### Bioconductor (18 paquetes clave)

| Paquete | Rol |
|---|---|
| DESeq2 | Expresión diferencial pseudobulk |
| clusterProfiler | Enriquecimiento GO |
| org.At.tair.db | Base de datos Arabidopsis |
| scater, scuttle | QC single-cell |
| SingleCellExperiment | Estructura de datos SC |
| SummarizedExperiment | Infraestructura Bioconductor |
| zellkonverter | Conversión R ↔ Python (h5ad) |
| BiocParallel | Paralelización en pipelines Bioc |

### GitHub (instalación directa)

| Paquete | Repositorio | Rol |
|---|---|---|
| SeuratDisk | mojaveazure/seurat-disk | Leer/escribir H5AD |
| DoubletFinder | chris-mcginnis-ucsf/DoubletFinder | Detección dobletes |
| SeuratWrappers | satijalab/seurat-wrappers | Integraciones externas |
| monocle3 | cole-trapnell-lab/monocle3 | Trayectorias celulares |
| Signac | — | Dependencia de SeuratWrappers |

### Legacy/Archivo

- `grr` (CRAN Archive) — compatibilidad con código antiguo

---

## CellRanger dentro del contenedor

- **Versión**: 9.0.1
- **Ruta**: `/opt/cellranger-9.0.1/`
- **Instalación**: Binario descargado y añadido al PATH durante el build
- **Rol**: Alineamiento de lecturas FASTQ → matriz de conteos (upstream del pipeline R)

---

## Opciones globales de R (Rprofile.site)

```r
options(repos = c(CRAN = "https://packagemanager.posit.co/cran/latest"))
options(Ncpus = parallel::detectCores())
```

---

## Backlinks

- [[ScRNASeq-Docker MOC]]
- [[ScRNASeq - Flujo de Datos]]
