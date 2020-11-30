#' @title Load and Process Drug Reference Profiles.
#' @description  This function allows user to load in the tissue specific drug rank matrix.
#' @details This function is a reverised version of get.cmap.ref from DrInsight package. The tissue specific drug rank matrix is tranformed from L1000data (GEO: GSE92742 and GSE70138) using PrepareReference function.
#' @param drug.response.path The local path and the name of the tissue specific drug rank matrix.
#' @param probe.to.genes A data.frame contains gene IDs (the IDs used in drug rank matrix) and official gene symbol. This files was automately generated with drug rank matrix.
#' @param drug.info A data.frame contains drug information. This file was automately generated with drug rank matrix.
#' @export


GetDrugRef = function(drug.response.path = NULL, probe.to.genes = NULL, drug.info = NULL){
  cat("\n")
  cat("\n")
  message("Loading CMap drug matrix. This may take some time ... \n")
  cmap.drug.rank = read.table(drug.response.path,row.names = 1, header = T, check.names = FALSE)
  cmap.drug.rank = cmap.drug.rank[probe.to.genes$ID,]
  rownames(cmap.drug.rank) = probe.to.genes$Gene.Symbol
  cmap.ref.profiles = list(drug.info = drug.info, drug.rank.matrix = cmap.drug.rank)
  return(cmap.ref.profiles)
}

