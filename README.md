# VOC-Targeted Breathomics Biomarker Pipeline

## 📖 Project Overview
This repository contains a comprehensive, end-to-end R pipeline for discovering and validating Volatile Organic Compound (VOC) biomarkers in exhaled breath. Originally designed for occupational health research (e.g., characterizing metabolic profiles in Coal Workers' Pneumoconiosis), the pipeline integrates high-dimensional metabolomics data with clinical covariates to identify stable diagnostic signatures.

### Pipeline Highlights
1. **Global Metabolic Profiling:** 3D Principal Component Analysis (PCA) and hierarchical clustering heatmaps (`ComplexHeatmap`).
2. **Robust Feature Selection:** A two-step screening process utilizing covariate-adjusted `limma` (False Discovery Rate control) followed by LASSO penalized regression (`glmnet`).
3. **Machine Learning Diagnostics:** Training and internal validation of Support Vector Machine (SVM), Random Forest (RF), and Firth Logistic Regression models, evaluated via ROC/AUC analysis.
4. **Clinical Multi-omics Integration:** Spearman correlation mapping and Cytoscape-like network visualizations (`ggraph`) linking selected VOCs with ELISA inflammatory markers (e.g., TNF-α, IL-6) and lung function indices (e.g., FVC, FEV1).

---

## 📂 Repository Structure: Code & Data

### 📊 Data Requirements
To run this pipeline, you need the following datasets formatted as `.csv` files:
* `PN_DE_data.csv`: Master dataset containing Sample IDs, Group labels (e.g., HC, DE, PN), clinical covariates (Age, BMI, Smoke_status, Complications), and raw/normalized VOC intensities.
* `VOC_preprocessed_R_output.csv`: Intermediate VOC matrix used for correlation analyses.
* `elisa_voc.csv`: Clinical dataset containing ELISA inflammatory marker concentrations.

*(Note: Due to patient privacy and data protection regulations, raw clinical data is not included in this public repository. Please use the provided mock templates or your own datasets following the exact column structures).*

### 💻 Core R Scripts
* `1_VOC_Biomarker_Discovery_Pipeline.R`: The main machine learning pipeline for subsetting groups (e.g., PN vs. HC), data normalization, Limma-LASSO feature selection, modeling, and plotting ROCs/Boxplots.
* `2_Integrated_Correlation_Analysis.R`: The clinical correlation script. It securely merges VOC data with ELISA and lung function data, handling encoding issues, computing FDR-adjusted Spearman correlations, and generating heatmaps and network graphs.

---

## ⚙️ Environment Configuration

This pipeline was developed and tested in **R (>= 4.1.0)**. 

### Required Packages
You will need several standard CRAN packages and Bioconductor tools. Run the following snippet in your R console to configure your environment:

```R
# Install CRAN packages
cran_pkgs <- c("tidyverse", "caret", "glmnet", "e1071", "randomForest", 
               "logistf", "pROC", "plotly", "ggpubr", "patchwork", 
               "writexl", "igraph", "ggraph", "tidygraph", "Hmisc", "reshape2")
new_cran <- cran_pkgs[!(cran_pkgs %in% rownames(installed.packages()))]
if(length(new_cran)>0) install.packages(new_cran)

# Install Bioconductor packages (for limma and ComplexHeatmap)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
bioc_pkgs <- c("limma", "ComplexHeatmap")
new_bioc <- bioc_pkgs[!(bioc_pkgs %in% rownames(installed.packages()))]
if(length(new_bioc)>0) BiocManager::install(new_bioc)
