#' @title Single Drug Selection.
#' @description  Select single drugs for every cell population by FDR and drug type, and summarize cell coverage for selected drugs.
#' @details Input raw drug repurosing result and return the top drugs with summary of cell coverage.
#' @param SC.integrated Single cell alignment result from SCalignment function.
#' @param Drug.data Drug repurosing result from drug.ident function.
#' @param Drug.FDR The FDR threshold to select drug. The default value is 0.1.
#' @param FDA.drug.only logical; if TRUE, will only return FDA-approved drugs.
#' @param Case An vector of case (diseased) samples.Only case sammples are involved in the calculation of coverage.
#' @return A data frame of selected drugs with summary of cell coverage.
#' @export
#' @import Seurat

TopDrug <- function(SC.integrated = SC.data,
                    Drug.data = Drug.ident.res,
                    Drug.FDR = 0.1,
                    FDA.drug.only = TRUE,
                    Case = NULL){
  data <- SC.integrated
  ##Cell proportion
  cells <- SC.integrated@meta.data
  if(length(Case)>0){
  cells <- subset(cells,sample %in% Case)
  }
  cells <- cells$celltype
  cell.count <- table(cells)
  cell.count <- cell.count[which(cell.count>3)]
  cells.freq <- round(100*cell.count/length(cells),2)

  ##Load drug data
  Drug.list <- data.frame()
  for(i in names(Drug.data)){
    Cd <- Drug.data[[i]]
    Cd <- subset(Cd, FDR<Drug.FDR)
    Drugs <- Cd$Drug.name
    if(FDA.drug.only==TRUE){
      Drugs <- intersect(Drugs,FDA.drug)
    }
    if(length(Drugs)>0){
      Cd <- subset(Cd, Drug.name %in% Drugs)
      temp <- data.frame(Drug=Cd$Drug.name,Cell.type=i,Cell.type.coverage=cells.freq[i],FDR=Cd$FDR,row.names = NULL)
      Drug.list <- rbind(Drug.list,temp)
    }
  }
  Drug.list <- Drug.list[order(Drug.list$FDR, decreasing = F),]
  Drug.list <- Drug.list[!duplicated(Drug.list),]
  Drug.coverage <- tapply(Drug.list$Cell.type.coverage, Drug.list$Drug,sum)
  temp.coverage <- Drug.coverage[Drug.list$Drug]
  Drug.list$Drug.coverage <- temp.coverage
  Drug.list <- Drug.list[,c(1:3,5,4)]
  Drug.list <- Drug.list[order(Drug.list$Drug.coverage, decreasing = T),]
  return(Drug.list)
}
