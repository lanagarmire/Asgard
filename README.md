## Asgard: A Single-cell Guided pipeline for Accurate Repurposing and combination of Drugs 
Using scRNA-seq data, Asgard repurposes mono-drugs for every single cell population and predicts personalized drug combinations to address the cellular heterogeneity of patients. 
## Installation
### Install dependency packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("SingleR","limma","cmapR")
install.packages('Seurat')
### Install Asgard
devtools::install_github("lanagarmire/Asgard")
