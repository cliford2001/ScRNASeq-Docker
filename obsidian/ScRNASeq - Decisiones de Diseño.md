---
tags: [scrnaseq, arquitectura, decisiones, quirks]
type: note
relacionado: [[ScRNASeq-Docker MOC]]
---

# ScRNASeq — Decisiones de Diseño

Justificaciones técnicas detrás de las elecciones del pipeline, incluyendo advertencias y comportamientos no obvios.

---

## Decisiones arquitectónicas

### 1. Separación funciones / orquestador

**Elección**: `ScRNA_Analysis_Functions.R` (biblioteca) separado de `scrnaseq_pipeline.R` (script).

**Por qué**: Permite reutilizar funciones en proyectos distintos, hacer pruebas independientes y personalizar workflows sin tocar el orquestador principal.

**Implicación práctica**: Para usar solo una función específica basta con:
```r
source("load_libraries.R")
source("ScRNA_Analysis_Functions.R")
```

---

### 2. CellBender como paso previo recomendado

**Elección**: El pipeline espera input en formato CellBender HDF5, no CellRanger directo.

**Por qué**: El RNA ambiente es especialmente problemático en tejidos vegetales (contaminación del buffer de protoplastación). CellBender remueve este ruido antes de que entre al pipeline.

**Fallback**: `USE_CELLBENDER = FALSE` en la configuración permite input desde CellRanger.

**Quirk técnico**: `load_cellbender_filtered_h5()` parsea manualmente los componentes CSR de la matriz sparse (data, indices, indptr, shape), en lugar de usar un lector HDF5 estándar. Más frágil pero más portable.

---

### 3. QC cloroplástico específico para plantas

**Elección**: Métricas separadas `percent.mt` (mitocondrial) y `percent.cp` (cloroplástico).

**Por qué**: Las células vegetales tienen cloroplastos funcionales; el RNA cloroplástico no indica células muertas como sí ocurre en animales.

**Para animales**: Pasar `cp_pattern=NULL` a `load_sample()` desactiva esta métrica.

**Patrones para Arabidopsis**:
- Mitocondrial: `^ATMG`
- Cloroplástico: `^ATCG`

---

### 4. Harmony en lugar de RPCA

**Elección**: `RunHarmony("orig.ident")` para corrección de batch.

**Por qué**: Opera en espacio PCA sin requerir features comunes entre muestras. Escala bien a 5+ muestras.

**Cuándo reconsiderar**: Efectos de batch severos o integración cross-especie → considerar scVI via SeuratWrappers.

**Consecuencia crítica**: Todos los steps posteriores usan `reduction = "harmony"`, no PCA.

---

### 5. Pseudo-replicación para DESeq2

**Elección**: `asignar_pseudoreplicados()` asigna células al azar a `n_reps=3` grupos.

**Por qué**: DESeq2 requiere réplicas biológicas. Experimentos de una sola réplica no pueden usarse directamente. Los pseudo-replicados son un compromiso aceptado en la literatura.

**Advertencia**: Infla estadísticamente los grados de libertad. Los resultados deben interpretarse con cautela.

---

### 6. Checkpoints RDS intermedios

**Elección**: Guardar objetos Seurat en puntos clave del pipeline.

**Por qué**: El pipeline tarda potencialmente días en datasets grandes. Los checkpoints permiten re-entrar desde un estado guardado.

**Cómo restaurar**:
```r
# Descomentados en scrnaseq_pipeline.R por defecto
pbmc_harmony <- readRDS("dir_04/pbmc_harmony_postharmony.rds")
```

---

### 7. Configuración hardcodeada en el script

**Elección**: No hay archivo de configuración externo (YAML/JSON/INI).

**Por qué**: Simplicidad y reproducibilidad — todos los parámetros están "congelados" en el control de versiones junto al código.

**Desventaja**: Cambiar un parámetro requiere editar el script directamente. No hay CLI.

---

### 8. Sección 12 solo interactiva

**Elección**: La sección de curación manual no puede ejecutarse con `source()` completo.

**Por qué**: Requiere inspección humana de figuras intermedias para decidir qué células reasignar.

**Consecuencia**: El pipeline completo no puede correr de forma completamente automatizada.

---

### 9. Cacheo de `FindAllMarkers`

**Elección**: `find_markers()` cachea el resultado en TSV y lo lee en ejecuciones subsiguientes.

**Por qué**: `FindAllMarkers()` es computacionalmente costoso en datasets grandes (~horas).

**Cómo forzar recalcular**: `find_markers(obj, force=TRUE)`

---

## Quirks y advertencias

### Seurat 5: `JoinLayers()` antes de subsetear

En Seurat 5, los datasets integrados almacenan matrices por capa separada. Subsetear sin unificar falla silenciosamente:

```r
# ANTES de subsetear o hacer operaciones en SECTION 10:
obj <- JoinLayers(obj)
```

---

### Tasa de dobletes: 7.5% vs. estándar 10x

El `expected_doublet_rate=0.075` es mayor que el estándar de 10x Genomics (~0.8% per 1.000 células). Puede sobre-filtrar células reales en algunos datasets.

**Ajustar según**: densidad de carga, protocolo de captura, y resultados previos.

---

### Variable global `output_dir`

Las funciones de guardado dependen de una variable global. Fácil de olvidar al inicio de cada sección:

```r
output_dir <- dir_06   # ¡Hacer esto antes de save_pdf()!
```

---

### Atlas de referencia GSE273033

Usado en `annotate_by_reference()` para transferencia de etiquetas independiente. Específico de Arabidopsis root. Para otros organismos u órganos se necesita un atlas distinto.

---

### Organización de resultados por sección

`create_pipeline_dirs()` crea `dir_00` a `dir_12`. La convención es actualizar `output_dir` al inicio de cada sección para que las figuras se guarden en el directorio correcto.

---

### Semilla 1807

`set.seed(1807)` se usa en todo el pipeline. El número parece arbitrario (posiblemente una fecha o ID personal). Cambiarla produce resultados ligeramente distintos en pasos estocásticos (UMAP, clustering, pseudo-replicados).

---

## Backlinks

- [[ScRNASeq-Docker MOC]]
- [[ScRNASeq - Configuración]]
- [[ScRNASeq - Flujo de Datos]]
