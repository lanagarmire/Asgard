#' @title Single-cell Alignment.
#' @description  Single-cell pairwise correspondence, clustering and annotation.
#' @details In Asgard, every cell type in the diseased sample is paired to that in the normal (or control) sample, according to “anchor” genes that are consistently expressed between diseased and normal cells. It then clusters paired cells and identifies differentially expressed genes between paired diseased and normal cells for every clustered single cell population. Alternatively, in a cell-type-based procedure, Asgard annotates cell type for every single cell population and identifies differentially expressed genes between paired diseased and normal cells for every type of cell.
#' @param SC.list A list contains single-cell data.
#' @param CellCycle logical; if TRUE, will do cell-cycle regression. The default value is TURE.
#' @param anchor.features The number of "anchor" genes used in single-cell pairwise correspondence. The default value is 2000.
#' @param by.CellType logical; if TRUE, will annotate cell type. The default value is FALSE.
#' @return A Seurat object of aligned cells.
#' @export
#' @import Seurat SingleR celldex

SCalignment <- function(SC.list=SC.list,CellCycle=TRUE,anchor.features=2000,by.CellType=FALSE){
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
    SC.integrated@meta.data$seurat_clusters <- SC.integrated@meta.data$seurat_clusters + 1

    ##Cell Type Annotation
    if(by.CellType == TRUE){
     data <- as.matrix(SC.integrated@assays$RNA@data)
     hpca.se <- celldex::HumanPrimaryCellAtlasData
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
      SC.integrated@meta.data
      cells.temp <- data.frame(cell.id=row.names(sub.data@meta.data),celltype=best.cell)
      new.cells <- rbind(new.cells,cells.temp)
     }
     cell.meta <- SC.integrated@meta.data
     cell.id <- rownames(cell.meta)
     row.names(new.cells) <- new.cells[,1]
     new.cells <- new.cells[cell.id,]
     SC.integrated@meta.data$celltype <- new.cells$celltype
    }else{
     SC.integrated@meta.data$celltype <- paste0("C",SC.integrated@meta.data$seurat_clusters)
    }
    return(SC.integrated)
}

