FROM rocker/r-ver:4.5

# ── System dependencies ────────────────────────────────────────────────────────
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev libssl-dev libxml2-dev \
    libfontconfig1-dev libharfbuzz-dev libfribidi-dev \
    libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \
    libgdal-dev libgeos-dev libproj-dev libsqlite3-dev libudunits2-dev \
    libhdf5-dev \
    libv8-dev \
    libgit2-dev libssh2-1-dev cmake make git wget curl \
    python3 python3-pip python3-venv \
    libglpk-dev libfftw3-dev libgsl-dev \
    && rm -rf /var/lib/apt/lists/*

# ── Python virtual environment + packages ─────────────────────────────────────
RUN python3 -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir \
    scanpy \
    scFates \
    palantir \
    pandas \
    numpy \
    scipy \
    scikit-learn \
    matplotlib \
    seaborn \
    cellbender

# ── R global options ───────────────────────────────────────────────────────────
RUN echo 'options(repos = c(CRAN = "https://packagemanager.posit.co/cran/__linux__/noble/latest"))' \
    >> /usr/local/lib/R/etc/Rprofile.site && \
    echo 'options(Ncpus = parallel::detectCores())' \
    >> /usr/local/lib/R/etc/Rprofile.site

# ── BiocManager ────────────────────────────────────────────────────────────────
RUN R -e "install.packages('BiocManager')"

# ── CRAN packages ─────────────────────────────────────────────────────────────
RUN R -e "install.packages(c( \
    'Matrix', 'sp', 'Rcpp', \
    'ggplot2', 'patchwork', 'dplyr', 'tibble', 'knitr', 'kableExtra', \
    'tidyverse', 'tidyr', 'readr', 'purrr', 'forcats', 'lubridate', 'stringr', \
    'cowplot', 'gridExtra', \
    'UpSetR', 'eulerr', 'ggvenn', 'VennDiagram', 'futile.logger', \
    'clustree', 'ggraph', \
    'harmony', \
    'hdf5r', \
    'data.table', 'remotes', \
    'RColorBrewer', 'viridis', 'viridisLite', \
    'ggrepel', 'ggbeeswarm', 'ggridges', 'ggforce', \
    'plotly', 'htmlwidgets', 'shiny', 'miniUI', \
    'future', 'future.apply', 'parallelly', 'progressr', \
    'igraph', 'tidygraph', 'graphlayouts', \
    'lme4', 'nlme', 'MASS', 'boot', 'survival', \
    'reticulate', \
    'rmarkdown', 'xfun', 'htmltools', \
    'jsonlite', 'httr', 'xml2', \
    'devtools', 'KernSmooth', 'fields', 'ROCR', 'R.utils' \
))"

# ── Bioconductor packages ──────────────────────────────────────────────────────
RUN R -e "BiocManager::install(c( \
    'BiocGenerics', 'Biobase', 'S4Vectors', 'IRanges', \
    'GenomicRanges', 'SummarizedExperiment', 'SingleCellExperiment', \
    'MatrixGenerics', 'matrixStats', \
    'DelayedArray', 'SparseArray', 'S4Arrays', \
    'XVector', 'Seqinfo', \
    'scuttle', 'scater', \
    'BiocNeighbors', 'BiocParallel', 'BiocSingular', \
    'ScaledMatrix', 'beachmat', \
    'DESeq2', \
    'zellkonverter', \
    'basilisk', \
    'GenomeInfoDb', 'Rsamtools', \
    'BiocManager' \
), ask = FALSE, update = FALSE)"

# ── Seurat ecosystem ───────────────────────────────────────────────────────────
RUN R -e "install.packages(c('Seurat', 'SeuratObject'))"

# ── GitHub: SeuratDisk ────────────────────────────────────────────────────────
RUN R -e "remotes::install_github('mojaveazure/seurat-disk')"

# ── GitHub: DoubletFinder (no está en CRAN) ───────────────────────────────────
RUN R -e "remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')"

# ── GitHub: SeuratWrappers ────────────────────────────────────────────────────
ARG GITHUB_PAT
RUN R -e "install.packages(c('Signac'))" && \
    R -e "Sys.setenv(GITHUB_PAT='${GITHUB_PAT}'); remotes::install_github('satijalab/seurat-wrappers', upgrade='never')"

# ── grr (archivado CRAN) + monocle3 ──────────────────────────────────────────
RUN R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/grr/grr_0.9.5.tar.gz', repos=NULL, type='source')"
RUN R -e "remotes::install_github('cole-trapnell-lab/monocle3')"

# ── CellRanger 10.0.0 ─────────────────────────────────────────────────────────
COPY cellranger-9.0.1.tar.gz /opt/
RUN tar -xzf /opt/cellranger-9.0.1.tar.gz -C /opt/ && \
    rm /opt/cellranger-9.0.1.tar.gz
ENV PATH="/opt/cellranger-9.0.1:$PATH"

WORKDIR /workspace

CMD ["R"]
