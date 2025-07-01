# install.R – Heraclomics dependencies installer
# Run with: source("install.R")

# ---- SETTINGS ----
options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# ---- CRAN PACKAGES ----
cran_packages <- c(
  "shiny", "shinyWidgets", "shinyBS", "shinyjs", "shinydashboard", "shinycssloaders", "shinyalert",
  "plotly", "DT", "ggplot2", "ggrepel", "cowplot", "viridis", "reshape2", "pheatmap",
  "dplyr", "tibble", "tidyverse", "logger", "parallel", "future", "future.apply", "doParallel",
  "Matrix", "openxlsx", "uwot", "NMF", "RColorBrewer", "htmlwidgets", "visNetwork"
)

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}
invisible(lapply(cran_packages, install_if_missing))

# ---- BIOCONDUCTOR PACKAGES ----
bioc_packages <- c(
  "Seurat", "SeuratWrappers", "SingleR", "celldex", "edgeR", "DESeq2",
  "clusterProfiler", "org.Hs.eg.db", "AUCell", "RcisTarget", "GENIE3", "RTN", "HGNChelper", "enrichplot", "msigdbr"
)

invisible(lapply(bioc_packages, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg)
}))

# ---- GITHUB PACKAGES ----
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")

# DoubletFinder
if (!requireNamespace("DoubletFinder", quietly = TRUE)) {
  remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")
}

# CellChat
if (!requireNamespace("CellChat", quietly = TRUE)) {
  remotes::install_github("sqjin/CellChat")
}

# monocle3
if (!requireNamespace("monocle3", quietly = TRUE)) {
  if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
  devtools::install_github("cole-trapnell-lab/monocle3")
}

# MAGMA.Celltyping
if (!requireNamespace("MAGMA.Celltyping", quietly = TRUE)) {
  remotes::install_github("nathan-geffen/MAGMA_Celltyping")
}

# MEGENA
if (!requireNamespace("MEGENA", quietly = TRUE)) {
  remotes::install_github("bcbio/MEGENA")
}

# ---- TENSORFLOW & KERAS ----
if (!requireNamespace("keras", quietly = TRUE)) install.packages("keras")
if (!requireNamespace("tensorflow", quietly = TRUE)) install.packages("tensorflow")
library(tensorflow)
try(tensorflow::install_tensorflow(), silent = TRUE)

# ---- OPTIONAL: Download cisTarget Databases ZIP from Google Drive ----
zip_url <- "https://drive.google.com/uc?export=download&id=1sUVVkgl_GFvfuxgQa4N81EvpdXn9x914"
zip_dest <- "Scripts/cisTarget_databases.zip"
unzip_dir <- "Scripts"

if (!dir.exists(file.path(unzip_dir, "cisTarget_databases"))) {
  message("📥 Downloading cisTarget_databases.zip from Google Drive...")
  download.file(zip_url, destfile = zip_dest, mode = "wb")
  unzip(zip_dest, exdir = unzip_dir)
  unlink(zip_dest)
  message("✅ Extracted cisTarget_databases to 'Scripts/'")
} else {
  message("✔️ cisTarget_databases folder already exists.")
}

# ---- DONE ----
cat("\n✅ All required packages have been installed.\n")
