 # 🧬 Mouse TF-SEM Single Cell Analysis

This repository contains **R** and **Python** code used to generate the figures in:

> *"Transgene-Free Generation of Post-Gastrulation Mouse Whole Embryo Models Derived Solely from Naïve ESCs and iPSCs"*,  
> *Cell Stem Cell, 2025*

The R code was executed using **R 4.3.1**.  
The Python code environment is described in [`scvi_env.yml`](scANVI/scvi_env.yml).

---

## 📂 R Scripts

### Basic Seurat Analyses
- `A_2iLIF_seurat_analysis.R` — 2iLIF cell culture sample  
- `B_24h_seurat_analysis.R` — Day 1 (24h) cell culture sample  
- `C_48h_seurat_analysis.R` — Day 2 (48h) cell culture sample  
- `D_96h_seurat_analysis.R` — Day 4 (96h) cell culture sample  
- `E75_seurat_analysis.R` — E7.5 samples  
- `E85_seurat_analysis.R` — E8.5 samples  

### Advanced Analyses
- `cell_population_analysis.R` —  
  Comprehensive analysis across all time points, including:
  - UMAP plots  
  - Cluster and cell annotation  
  - Dotplot expression  
  - Feature expression  
  - Triangle analysis (2D projection of 3D data)  

- `cell_populations_on_ref_atlas.R` —  
  Projection of all time points onto the Proks *et al.* reference atlas.

- `monocole3_analysis.R` —  
  Pseudotime analysis of the 48h (Day 2) sample using Monocle3.

### Data Availability
rds files can be found in ncbi GEO:  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE294766

---

## 🧠 scANVI (Python)

Located in the `scANVI/` folder:
```
scANVI/
├── scvi_env.yml # Conda environment file
├── config_paper_3L_35_128.in # Configuration for training
├── config_SEM85_test_query.in # Configuration for transfer learning
├── learn_with_scVI.py # scVI/scANVI training script
├── test_query_with_model.py # Transfer learning on query data
└── utils.py # Utility functions
```
### Data Availability
Train data can be downloaded from zenodo: 
https://zenodo.org/records/15789211?token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6IjE1N2YxYmQ2LTViZDctNDllMy1hMmEzLWM1NTM4OGMyNjIxYyIsImRhdGEiOnt9LCJyYW5kb20iOiIyNGM2OTgwMTBmNmUwYzBkZDY0ZmRmNDUyNGJiZDUxZSJ9.daS_o360CcS3vguS06Ck2TZl1QfqBa-svmtxlV5nXSlso_R5zE57J4wjwYJBt8vyZl4uIHsRazliz51dBYzGQg
---

## 📌 Notes

- Ensure all required packages are installed (use `scvi_env.yml` for Python).
- Reference atlas data should be downloaded separately.
- If you use this code, please cite the original paper.

---
