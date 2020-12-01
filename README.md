## Asgard: A Single-cell Guided pipeline for Accurate Repurposing of Drugs 
Using scRNA-seq data, Asgard repurposes mono-drugs for every single cell population and predicts personalized drug combinations to address the cellular heterogeneity of patients. 
## Installation
### Install dependency packages
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    
BiocManager::install(c("SingleR","limma","cmapR"))

install.packages('Seurat')
```
### Install Asgard
```
devtools::install_github("lanagarmire/Asgard")
```
## Prepare Drug Referecne
### Step 1
#### Download L1000 Connectivity Map perturbational profiles GSE70138 and GSE92742 from GEO

[GSE70138_Broad_LINCS_cell_info_2017-04-28.txt](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/GSE70138_Broad_LINCS_cell_info_2017-04-28.txt.gz)

[GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx.gz)

[GSE70138_Broad_LINCS_sig_info_2017-03-06.txt](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/GSE70138_Broad_LINCS_sig_info_2017-03-06.txt.gz)

[GSE70138_Broad_LINCS_gene_info_2017-03-06.txt](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/GSE70138_Broad_LINCS_gene_info_2017-03-06.txt.gz)

[GSE92742_Broad_LINCS_cell_info.txt](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_cell_info.txt.gz)

[GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx.gz)

[GSE92742_Broad_LINCS_sig_info.txt](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_sig_info.txt.gz)

### Step 2 
#### Generate drug reference panel for lung from GSE70138 and GSE92742
Unzip downloaded files and run the following code:
```
library('Asgard')

PrepareReference(cell.info="GSE70138_Broad_LINCS_cell_info_2017-04-28.txt",
                 gene.info="GSE70138_Broad_LINCS_gene_info_2017-03-06.txt",
                 GSE70138.sig.info = "GSE70138_Broad_LINCS_sig_info_2017-03-06.txt",
                 GSE92742.sig.info = "GSE92742_Broad_LINCS_sig_info.txt",
                 GSE70138.gctx = "GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx",
                 GSE92742.gctx = "GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx",
                 Output.Dir = "DrugReference/"
)

```
Please use '?PrepareReference' for more help.

