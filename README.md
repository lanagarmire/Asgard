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
### Load Asgard
```
library('Asgard')
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
#### Generate tissue specific drug references from GSE70138 and GSE92742
Unzip downloaded files, revise the Your_local_path and run the following code:
```
library('Asgard')

PrepareReference(cell.info="Your_local_path/GSE70138_Broad_LINCS_cell_info_2017-04-28.txt",
                 gene.info="Your_local_path/GSE70138_Broad_LINCS_gene_info_2017-03-06.txt",
                 GSE70138.sig.info = "Your_local_path/GSE70138_Broad_LINCS_sig_info_2017-03-06.txt",
                 GSE92742.sig.info = "Your_local_path/GSE92742_Broad_LINCS_sig_info.txt",
                 GSE70138.gctx = "Your_local_path/GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx",
                 GSE92742.gctx = "Your_local_path/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx",
                 Output.Dir = "Your_local_path/DrugReference/"
)

```
Please use '?PrepareReference' for more help.

## Drug Repurposing
### Step 1
#### Load single-cell RNA-seq data
```
library('Seurat')
#Case sample
data<-read.table(file="Your_local_path/Case_Expression_Matrix.txt",header = T,check.names=FALSE)
Case <- CreateSeuratObject(counts = data, project = "Demo", min.cells = 3, min.features = 200,meta.data=data.frame(cell=colnames(data),sample="Case"))
#Control sample
data<-read.table(file="Your_local_path/Control_Expression_Matrix.txt",header = T,check.names=FALSE)
Control <- CreateSeuratObject(counts = data, project = "Demo", min.cells = 3, min.features = 200,meta.data=data.frame(cell=colnames(data),sample="Control"))

```

#### Single-cell alignment
```
library('Asgard')
#Creat data list
SC.list<-list(Case=Case,Control=Control)
#Single-cell alignment without cell type annotation
SC.data<-SCalignment(SC.list,CellCycle=TRUE,anchor.features=2000,by.CellType=FALSE)
#Single-cell alignment with cell type annotation (optional)
SC.data<-SCalignment(SC.list,CellCycle=TRUE,anchor.features=2000,by.CellType=TRUE)
#Visualize alignment result
DimPlot(SC.data, reduction = "umap", split.by = "sample", label = TRUE)
```

#### Single-cell comparison
```
#Case sample names
Case.samples=c("Case")
#Control sample names
Control.samples=c("Control")
#Get differential gene expression profiles for every cell type (or cluster if without annotation)
Gene.list<-GetGene(SC.integrated=SC.data,Case=Case.samples,Control=Control.samples,min.cells=3)

```

#### Mono-drug repurposing for every cell type
```
#Load tissue specific drug reference
my_gene_info<-read.table(file="Your_local_path/DrugReference/breast_gene_info.txt",sep="\t",header = T,quote = "")
my_drug_info<-read.table(file="Your_local_path/DrugReference/breast_drug_info.txt",sep="\t",header = T,quote = "")
drug.ref.profiles = GetDrugRef(drug.response.path = 'Your_local_path/DrugReference/breast_rankMatrix.txt',
                                 probe.to.genes = my_gene_info, drug.info = my_drug_info)
#Repurpose mono-drugs for every cell type                               
Drug.ident.res = GetDrug(gene.data = Gene.list, drug.ref.profiles = drug.ref.profiles, repurposing.unit = "drug", connectivity = "negative", drug.type = "FDA")
```
Use ?GetDrug for more help

#### Drug combination for every cell type
```
GSE92742.gctx.path="Your_local_path/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx"
GSE70138.gctx.path="Your_local_path/GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx"
Drug.combinations<-DrugCombination(SC.integrated=SC.data,
                      Gene.data=Gene.list,
                      Drug.data=Drug.ident.res,
                      Drug.FDR=0.1,
                      FDA.drug.only=TRUE,
                      Combined.drugs=2,
                      Case=Case.samples,
                      Tissue="breast",
                      GSE92742.gctx=GSE92742.gctx.path,
                      GSE70138.gctx=GSE70138.gctx.path)
```

#### Select mono-drug therapies
```
Final.drugs<-TopDrug(SC.integrated=SC.data,
                   Drug.data=Drug.ident.res,
                   Drug.FDR=0.1,
                   FDA.drug.only=TRUE,
                   Case=Case.samples
)
```

#### Select drug combination therapies
```
Final.combinations<-TopCombination(Drug.combination=Drug.combinations,
                   Combination.FDR=0.1,
                   Min.combination.score=1
)
```
