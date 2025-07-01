# 🧬 Heraclomics

Heraclomics is an interactive, no-code R Shiny application designed to streamline single-cell RNA-seq data analysis for researchers without programming expertise. It guides users from raw data preprocessing to advanced analyses like gene regulatory networks and cell-cell communication.

---

## 🚀 Key Features

- **End-to-end workflow**: Quality control, clustering, annotation, and downstream analysis  
- **Input flexibility**: Supports 10X Genomics, Seurat objects, or raw count matrices  
- **Advanced modules**:  
  - Doublet detection & batch correction  
  - Annotation (SingleR/scType)  
  - Differential expression, SCENIC, CellChat, Monocle3  
- **Interactive outputs**: Plotly visuals, downloadable tables/figures  

## ⚙️ System Requirements

- R ≥ 4.0.0  
- Bioconductor ≥ 3.12  
- 8GB+ RAM recommended for large datasets  

---

## 🛠️ Setup

Run the provided installation script to install all dependencies:

```r
source("install.R")
```

You can also install required packages manually from CRAN, Bioconductor, or GitHub.

---

## ▶️ Launching the App

```r
shiny::runApp("app.R")
```

Or via terminal:

```bash
Rscript app.R
```

---

## 🧬 RcisTarget Database Setup

> These databases are not installed automatically when download from Github due to their size.
You can dowload and create below structure

### 📁 Directory structure

```
Scripts/cisTarget_databases/
├── hgnc/
│   ├── hg19-500bp.........
│   ├── hg19-10kbp.........
│   └── 
├── mgi/
│   └── ...
└── dmel/
    └── ...
```

Download from: [https://resources.aertslab.org/cistarget/]

Or use the automatic download(*Recomended) built into `install.R`.

---



---


