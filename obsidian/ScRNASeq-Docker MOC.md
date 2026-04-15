---
tags: [scrnaseq, bioinformatics, docker, MOC]
type: MOC
---

# ScRNASeq-Docker — Mapa de Contenido

Pipeline reproducible de análisis de single-cell RNA-seq, contenedorizado con Docker. Desarrollado originalmente para *Arabidopsis thaliana* bajo estrés de nitrógeno, adaptable a otros organismos.

> **Repositorio**: https://github.com/cliford2001/ScRNASeq-Docker  
> **Directorio local**: `/Users/matias/ScRNASeq-Docker`

---

## Notas del proyecto

- [[ScRNASeq - Setup y Docker]] — Dockerfile, docker-compose, instalación, cómo correr el contenedor
- [[ScRNASeq - Flujo de Datos]] — De FASTQ a resultados: cada sección del pipeline paso a paso
- [[ScRNASeq - Funciones R]] — Referencia de las 50+ funciones en `ScRNA_Analysis_Functions.R`
- [[ScRNASeq - Configuración]] — Parámetros clave, cómo adaptar a otro organismo
- [[ScRNASeq - Decisiones de Diseño]] — Por qué se hicieron las elecciones técnicas, quirks y advertencias

---

## Resumen rápido

| Aspecto | Detalle |
|---|---|
| Lenguaje principal | R 4.5 (Seurat 5.4.0) |
| Lenguaje secundario | Python 3 (scanpy, palantir, scFates) |
| Input | HDF5 de CellBender o CellRanger |
| Output | RDS, CSVs, PDFs, `.h5ad` para scanpy |
| Organismo por defecto | *Arabidopsis thaliana* |
| Secciones del pipeline | 12 (SECTION 0–12) |
| Funciones de biblioteca | ~50 funciones / 2.764 líneas |
| Imagen Docker base | `rocker/r-ver:4.5` |
| Reproducibilidad | Contenedor + `set.seed(1807)` + paquetes versionados |

---

## Archivos clave del repo

```
ScRNASeq-Docker/
├── Dockerfile                      # Imagen del entorno
├── docker-compose.yml              # Orquestación del contenedor
├── load_libraries.R                # Carga de paquetes R
├── ScRNA_Analysis_Functions.R      # Biblioteca de funciones (~2.764 líneas)
├── custom_seurat.R                 # Visualizaciones Seurat personalizadas
├── scrnaseq_pipeline.R             # Orquestador end-to-end (~896 líneas)
├── arabidopsis_scrna.Rmd           # Protocolo documentado
└── README.md                       # Referencia de funciones
```

---

## Algoritmos clave

- **Harmony** → corrección de efectos de batch
- **DoubletFinder** → detección de dobletes
- **DESeq2** → expresión diferencial (pseudobulk)
- **Leiden/Louvain** → clustering
- **clusterProfiler / enrichGO** → enriquecimiento GO
- **CellBender** → remoción de RNA ambiente (pre-pipeline)
- **CellRanger 9.0.1** → procesamiento de lecturas (pre-pipeline)
