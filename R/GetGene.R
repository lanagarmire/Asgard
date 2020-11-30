#' @title Get Differntially Expressed Genes.
#' @description  Identify differnential gene expression profiles for case cells compared to control cells.
#' @details This function identifies differentially expressed genes between paired case and control cells for every type of cell. If there is no cell type annotation, it uses every cluster as a cell type.
#' @param SC.integrated A Seurat object of aligned cells from SCalignment function.
#' @param Case A vector contains names of case samples.
#' @param Control A vector contains names of control samples.
#' @param min.cells The minimum number of cells for a cell type. A cell type is omitted if it has less cells than the minimum number.
#' @return A list of differnential gene expression profiles for every cell type.
#' @export
#' @import Seurat limma


GetGene <- function(SC.integrated=SC.integrated,Case=Case,Control=Control,min.cells=3){
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
    return(Gene.list)
}

