#' @title Sinlge-cell Plasticity.
#' @description  It determines the plasticity of each cell type.
#' @details This function estimate the entropy of every cell in the case samples. For each cell type, it use the median entropy value as the plasticity of each cell type.
#' @param SC.integrated A Seurat object of aligned single cells from SCalignment function.
#' @param Case A vector contains names of case samples.
#' @return A data frame of plasticity, normailized plasticity and cell type coverage.
#' @export

SCplasticity <- function (SC.integrated = SC.data, Case=NULL) 
{
  if(length(Case)>0){
    SC.integrated <- subset(SC.integrated, sample %in% Case)
  }else{
    SC.integrated <- SC.integrated
  }
  SC.meta <- SC.integrated@meta.data
  expr.data <- as.matrix(SC.integrated@assays$RNA@counts)
  
  #Entorpy-based Plasticity
  probs   <- t(t(expr.data)/apply(expr.data,2,sum))
  probs[is.na(probs)] <- 0
  log.probs <- log(probs)
  log.probs[which(is.infinite(log.probs))] <- 0
  SC.meta$cell.entropy <- -apply(probs*log.probs/log(nrow(expr.data)),2,sum)
  SC.entropy <- tapply(SC.meta$cell.entropy, SC.meta$celltype, median)
  SC.entropy <- data.frame(Cell.Type=row.names(SC.entropy),Plasticity=SC.entropy)
  rm(expr.data)
  rm(log.probs)
  rm(probs)
  
  #Normalize Plasticity
  SC.entropy$Normalized.Plasticity=(SC.entropy$Plasticity-min(SC.entropy$Plasticity))/(max(SC.entropy$Plasticity)-min(SC.entropy$Plasticity))
  
  #Population Size
  Cluster.cell.rate <- table(SC.meta$celltype)/nrow(SC.meta)
  SC.entropy$Coverage <- 100*Cluster.cell.rate[row.names(SC.entropy)]
  
  return(SC.entropy)
}

