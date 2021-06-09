## Asgard: A Single-cell Guided pipeline for Accurate Repurposing of Drugs 
Using scRNA-seq data, Asgard repurposes mono-drugs for every single cell population and predicts personalized drug combinations to address the cellular heterogeneity of patients. 
## System Requirements
### Hardware requirements
Asgard package requires only a standard computer with enough RAM (>64GB) to support the in-memory operations.
### Software requirements
#### OS requirements
The package has been tested on the following systems:
```
Windows 10
CentOS Linux 7
```
#### R package dependencies
```
Seurat
limma
cmapR
SingleR
celldex
```
## Installation
#### Install devtools if you don't have it
```
install.packages('devtools')
```
#### Install dependency packages
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    
BiocManager::install(c("SingleR","limma","cmapR","celldex"))

install.packages('Seurat')

#If you can't install a package with above commands, try to download the gz file and install it locally.

#Take celldex package as an example:

#Downlaod the source package of celldex in linux
wget https://bioconductor.org/packages/release/data/experiment/src/contrib/celldex_1.0.0.tar.gz

#Start R
R

#Install celldex from the local source package
install.packages('celldex_1.0.0.tar.gz')

#Note: some dependency packages require R version newer than 4.0

```
#### Install Asgard
```
devtools::install_github("lanagarmire/Asgard")
```
#### Load Asgard
```
library('Asgard')
```
## Prepare Drug Referecne Library
#### Step 1
#### Download L1000 Connectivity Map perturbational profiles GSE70138 and GSE92742 from GEO
<p>Method 1: click file names below </p>

[GSE70138_Broad_LINCS_cell_info_2017-04-28.txt](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/GSE70138_Broad_LINCS_cell_info_2017-04-28.txt.gz)

[GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx.gz)

[GSE70138_Broad_LINCS_sig_info_2017-03-06.txt](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/GSE70138_Broad_LINCS_sig_info_2017-03-06.txt.gz)

[GSE70138_Broad_LINCS_gene_info_2017-03-06.txt](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/GSE70138_Broad_LINCS_gene_info_2017-03-06.txt.gz)

[GSE92742_Broad_LINCS_cell_info.txt](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_cell_info.txt.gz)

[GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx.gz)

[GSE92742_Broad_LINCS_sig_info.txt](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_sig_info.txt.gz)

or Method 2: run following commands in linux
```
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/GSE70138_Broad_LINCS_cell_info_2017-04-28.txt.gz
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx.gz
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/GSE70138_Broad_LINCS_sig_info_2017-03-06.txt.gz
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/GSE70138_Broad_LINCS_gene_info_2017-03-06.txt.gz
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_cell_info.txt.gz
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx.gz
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_sig_info.txt.gz
```
#### Step 2 
#### Generate tissue specific drug references from GSE70138 and GSE92742
Unzip downloaded files, revise the Your_local_path and run the following code:
```
library('Asgard')

#Please replace Your_local_path with your real local folder

PrepareReference(cell.info="Your_local_path/GSE70138_Broad_LINCS_cell_info_2017-04-28.txt",
                 gene.info="Your_local_path/GSE70138_Broad_LINCS_gene_info_2017-03-06.txt",
                 GSE70138.sig.info = "Your_local_path/GSE70138_Broad_LINCS_sig_info_2017-03-06.txt",
                 GSE92742.sig.info = "Your_local_path/GSE92742_Broad_LINCS_sig_info.txt",
                 GSE70138.gctx = "Your_local_path/GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx",
                 GSE92742.gctx = "Your_local_path/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx",
                 Output.Dir = "Your_local_path/DrugReference/"
)

#Note: the file names here maybe different after unzipping.
#Please note that it takes more than one hour to produce drug references in a standard computer with RAM>64GB.
```
Please use '?PrepareReference' for more help.

## Drug Repurposing
#### Step 1
#### Load single-cell RNA-seq data
<p>Download datasets GSE113197 and GSE123926 from GEO before running this script.</p>

Human Breast Cancer Epithelial Cells (GSE123926): [GSE123926_RAW.tar](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE123926&format=file) 

Normal Human Breast Epithelial Cells (GSE113197): [GSE113197_RAW.tar](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE113197&format=file) 

```
library('Seurat')

#Load normal sample Ind5 from GSE113197 dataset
celltype<-read.table(file="Normal_celltype.txt",header = T,check.names=FALSE)
data<-read.table(file="GSM3099847_Ind5_Expression_Matrix.txt",header = T,check.names=FALSE)
row.names(data)<-data[,1]
data<-data[,-1]
celltype2<-subset(celltype,sample=="Ind5" & celltype %in% c("Luminal_L2_epithelial_cells","Luminal_L1.1_epithelial_cells", "Luminal_L1.2_epithelial_cells", "Basal_epithelial_cells"))
common <- intersect(colnames(data), rownames(celltype2))
data<-data[,common]
Epithelial2 <- CreateSeuratObject(counts = data, project = "Epithelial", min.cells = 3, min.features = 200,meta.data=data.frame(celltype2,cell=colnames(data),type="Normal"))

#Load normal sample Ind6 from GSE113197 dataset
data<-read.table(file="GSM3099848_Ind6_Expression_Matrix.txt",header = T,check.names=FALSE)
row.names(data)<-data[,1]
data<-data[,-1]
celltype3<-subset(celltype,sample=="Ind6" & celltype %in% c("Luminal_L2_epithelial_cells","Luminal_L1.1_epithelial_cells", "Luminal_L1.2_epithelial_cells", "Basal_epithelial_cells"))
common <- intersect(colnames(data), rownames(celltype3))
data<-data[,common]
Epithelial3 <- CreateSeuratObject(counts = data, project = "Epithelial", min.cells = 3, min.features = 200,meta.data=data.frame(celltype3,cell=colnames(data),type="Normal"))

#Load normal sample Ind7 from GSE113197 dataset
data<-read.table(file="GSM3099849_Ind7_Expression_Matrix.txt",header = T,check.names=FALSE)
row.names(data)<-data[,1]
data<-data[,-1]
celltype4<-subset(celltype,sample=="Ind7" & celltype %in% c("Luminal_L2_epithelial_cells","Luminal_L1.1_epithelial_cells", "Luminal_L1.2_epithelial_cells", "Basal_epithelial_cells"))
common <- intersect(colnames(data), rownames(celltype4))
data<-data[,common]
Epithelial4 <- CreateSeuratObject(counts = data, project = "Epithelial", min.cells = 3, min.features = 200,meta.data=data.frame(celltype4,cell=colnames(data),type="Normal"))

#Load cancer sample PDX110 from GSE123926 dataset
TNBC_PDX.data<- Read10X(data.dir = "GSM3516947_PDX110")
TNBC.PDX2 <- CreateSeuratObject(counts = TNBC_PDX.data, project = "TNBC", min.cells = 3, min.features = 200, meta.data=data.frame(row.names=colnames(TNBC_PDX.data), cell=colnames(TNBC_PDX.data), sample="PDX-110",type="TNBC.PDX"))

#Load cancer sample PDX322 from GSE123926 dataset
TNBC_PDX.data<- Read10X(data.dir = "GSM3516948_PDX322")
TNBC.PDX3 <- CreateSeuratObject(counts = TNBC_PDX.data, project = "TNBC", min.cells = 3, min.features = 200, meta.data=data.frame(row.names=colnames(TNBC_PDX.data), cell=colnames(TNBC_PDX.data), sample="PDX-332",type="TNBC.PDX"))


```

#### Step 2
#### Single-cell alignment
```
SC.list<-list(TNBC.PDX2=TNBC.PDX2,TNBC.PDX3=TNBC.PDX3,Epithelial2=Epithelial2,Epithelial3=Epithelial3,Epithelial4=Epithelial4)
CellCycle=TRUE #Set it TRUE if you want to do Cell Cycle Regression
anchor.features=2000

for (i in 1:length(SC.list)) {
     SC.list[[i]] <- NormalizeData(SC.list[[i]], verbose = FALSE)
     SC.list[[i]] <- FindVariableFeatures(SC.list[[i]], selection.method = "vst",
                           nfeatures = anchor.features, verbose = FALSE)
    }
    SC.anchors <- FindIntegrationAnchors(object.list = SC.list,anchor.features = anchor.features, dims = 1:15)
    SC.integrated <- IntegrateData(anchorset = SC.anchors, dims = 1:15)
    DefaultAssay(SC.integrated) <- "integrated"
    if(CellCycle){
    ##Cell Cycle Regression
    s.genes <- cc.genes$s.genes
    g2m.genes <- cc.genes$g2m.genes
    SC.integrated <- CellCycleScoring(SC.integrated, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
    SC.integrated <- ScaleData(SC.integrated, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(SC.integrated))
    SC.integrated <- RunPCA(SC.integrated, npcs = 15, verbose = FALSE)
    }else{
     ##Run the standard workflow for visualization and clustering
     SC.integrated <- ScaleData(SC.integrated, verbose = FALSE)
     SC.integrated <- RunPCA(SC.integrated, npcs = 15, verbose = FALSE)
    }
    ##t-SNE and Clustering
    SC.integrated <- RunUMAP(SC.integrated, reduction = "pca", dims = 1:15)
    SC.integrated <- FindNeighbors(SC.integrated, reduction = "pca", dims = 1:15)
    SC.integrated <- FindClusters(SC.integrated, algorithm = 1, resolution = 0.4)

    ##Cell Type Annotation, set by.CellType=TRUE if you want to annotate cell  type.
    by.CellType=FALSE
    if(by.CellType == TRUE){
     data <- as.matrix(SC.integrated@assays$RNA@data)
     hpca.se <- HumanPrimaryCellAtlasData()
     pred.hpca <- SingleR(test = data, ref = hpca.se, assay.type.test=1, labels = hpca.se$label.main)
     cell.label <- data.frame(row.names = row.names(pred.hpca),celltype=pred.hpca$labels)
     if(length(SC.integrated@meta.data$celltype)>0){
      SC.integrated@meta.data$celltype <- cell.label$celltype
     }else{
       SC.integrated@meta.data <- cbind(SC.integrated@meta.data,cell.label)
     }
     new.cells <- data.frame()
     for(i in unique(SC.integrated$seurat_clusters)){
      sub.data <- subset(SC.integrated,seurat_clusters==i)
      temp <- table(sub.data@meta.data$celltype)
      best.cell <- names(which(temp==temp[which.max(temp)]))
      cells.temp <- data.frame(cell.id=row.names(sub.data@meta.data),celltype=best.cell)
      new.cells <- rbind(new.cells,cells.temp)
     }
     cell.meta <- SC.integrated@meta.data
     cell.id <- rownames(cell.meta)
     row.names(new.cells) <- new.cells[,1]
     new.cells <- new.cells[cell.id,]
     SC.integrated@meta.data$celltype <- new.cells$celltype
    }else{
     SC.integrated@meta.data$celltype <- paste0("C",as.numeric(SC.integrated@meta.data$seurat_clusters))
    }

#Change sample names
sample<-SC.data@meta.data$sample
sample[which(sample=="Ind5")]<-"Normal1"
sample[which(sample=="Ind6")]<-"Normal2"
sample[which(sample=="Ind7")]<-"Normal3"
SC.data@meta.data$sample<-sample

#Visualize alignment result
DimPlot(SC.data, reduction = "umap", split.by = "sample",group.by = "celltype")
```
#### Step 3
#### Single-cell comparison
```
#Case sample names
Case=c("PDX-110","PDX-332")

#Control sample names
Control=c("Normal1","Normal2","Normal3")


#Get differential gene expression profiles for every cell type (or cluster if without annotation) from Limma
library('limma')
DefaultAssay(SC.integrated) <- "RNA"
set.seed(123456)
Gene.list <- list()
C_names <- NULL
for(i in unique(SC.integrated@meta.data$celltype)){
     Idents(SC.integrated) <- "celltype"
     c_cells <- subset(SC.integrated, celltype == i)
     Idents(c_cells) <- "type"
     Samples=c_cells@meta.data
     Controlsample <- row.names(subset(Samples,sample %in% Control))
     Casesample <- row.names(subset(Samples,sample %in% Case))
     if(length(Controlsample)>min.cells & length(Casesample)>min.cells){
      expr <- as.matrix(c_cells@assays$RNA@data)
      new_expr <- as.matrix(expr[,c(Casesample,Controlsample)])
      new_sample <- data.frame(Samples=c(Casesample,Controlsample),type=c(rep("Case",length(Casesample)),rep("Control",length(Controlsample))))
      row.names(new_sample) <- paste(new_sample$Samples,row.names(new_sample),sep="_")
      expr <- new_expr
      bad <- which(rowSums(expr>0)<3)
      expr <- expr[-bad,]
      mm <- model.matrix(~0 + type, data = new_sample)
      fit <- lmFit(expr, mm)
      contr <- makeContrasts(typeCase - typeControl, levels = colnames(coef(fit)))
      tmp <- contrasts.fit(fit, contrasts = contr)
      tmp <- eBayes(tmp)
      C_data <- topTable(tmp, sort.by = "P",n = nrow(tmp))
      Gene.list[[i]] <- C_data
      C_names <- c(C_names,i)
      C_data_for_drug <- data.frame(geneSymbol=row.names(C_data),score=C_data$t)
     }
}
names(Gene.list) <- C_names

```

#### Step 4
#### Mono-drug repurposing for every cell type
```
library('Asgard')

#Please replace Your_local_path with your real local folder

#Load tissue specific drug reference produced by PrepareReference function as mentioned above
my_gene_info<-read.table(file="Your_local_path/DrugReference/breast_gene_info.txt",sep="\t",header = T,quote = "")
my_drug_info<-read.table(file="Your_local_path/DrugReference/breast_drug_info.txt",sep="\t",header = T,quote = "")
drug.ref.profiles = GetDrugRef(drug.response.path = 'Your_local_path/DrugReference/breast_rankMatrix.txt',
                               probe.to.genes = my_gene_info, 
                               drug.info = my_drug_info)

#Repurpose mono-drugs for every cell type                               
Drug.ident.res = GetDrug(gene.data = Gene.list, 
                        drug.ref.profiles = drug.ref.profiles, 
                        repurposing.unit = "drug", 
                        connectivity = "negative", 
                        drug.type = "FDA")
                       
```
Use '?GetDrug' for more help
#### Step 5
#### Drug combination analysis
```
#Please replace Your_local_path with your real local folder

GSE92742.gctx.path="Your_local_path/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx"
GSE70138.gctx.path="Your_local_path/GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx"
Drug.combinations<-DrugCombination(SC.integrated=SC.data,
                      Gene.data=Gene.list,
                      Drug.data=Drug.ident.res,
                      Drug.FDR=0.1,
                      FDA.drug.only=TRUE,
                      Combined.drugs=2,
                      Case=Case,
                      Tissue="breast",
                      GSE92742.gctx=GSE92742.gctx.path,
                      GSE70138.gctx=GSE70138.gctx.path)
```
Please use '?DrugCombination' for more help.
#### Step 6
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
Demo codes using real datasets are available at: https://github.com/lanagarmire/Single-cell-drug-repositioning

If you have further questions or comments, please contact Dr.Bing He: hbing@med.umich.edu
