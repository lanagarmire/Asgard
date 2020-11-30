#' @title Mono-drug Repurposing.
#' @description  Identify mono-drug therapy for every cell type.
#' @details This function allows user to use the differential expression data of every case cell type to query against reference drug response profiles.This function is a reverised version of drug.identification from DrInsight package.
#' @param drug.ref.profiles A list contains tissue specific drug reference Profiles from GetDrugRef function.
#' @param repurposing.unit The parameter of either "treatment" or "drug", which indicates if user want the function to test drug repurposing p value at treatment level or drug level. The default is "treatment", which treats the drug data from different cell lines separately.
#' @param CEG.threshold The p value threshold to select the consistently differential expressed genes (CEGs). The default value is 0.05.
#' @param connectivity The type of connectivity, either "negative" or "positive". Negative connectivity is used when the query data is the differential scores from disease data, and it will repurpose drugs that can potentially reverse the query disease phenotype. Positive connectivity is used when the query data is from a drug profile, and it will return the drugs that are similar to the query drug. The default value is "negative".
#' @param drug.type The parameter of either "FDA" or "compounds" or "all", which indicates if user want the function to identify FDA-approved drugs or compounds or both, respectively.The default value is "FDA".
#' @return A list of mono-drugs for every cell type.
#' @export


GetDrug = function(gene.data = NULL,
                    drug.ref.profiles = NULL,
                    repurposing.unit = "drug",
                    CEG.threshold = 0.05,
                    connectivity = "negative",
                    drug.type="FDA"){
  if(drug.type=="FDA"){
    Drug.info <- drug.ref.profiles$drug.info
    Drug.info$temp_name <- gsub("_.*","",Drug.info$cmap_name)
    Drug.info <- subset(Drug.info, temp_name %in% FDA.drug)
    Drug.info <- Drug.info[,colnames(drug.ref.profiles$drug.info)]
    drug.ref.profiles$drug.rank.matrix <- drug.ref.profiles$drug.rank.matrix[,Drug.info$instance_id]
    drug.ref.profiles$drug.info <- Drug.info
  }else if(drug.type=="compounds"){
    Drug.info <- drug.ref.profiles$drug.info
    Drug.info$temp_name <- gsub("_.*","",Drug.info$cmap_name)
    Drug.info <- subset(Drug.info, !(temp_name %in% FDA.drug))
    Drug.info <- Drug.info[,colnames(drug.ref.profiles$drug.info)]
    drug.ref.profiles$drug.rank.matrix <- drug.ref.profiles$drug.rank.matrix[,Drug.info$instance_id]
    drug.ref.profiles$drug.info <- Drug.info
  }
  res.list <- list()
  for(ci in 1:length(names(gene.data))){
      query.data <- data.frame(geneSymbol=row.names(gene.data[[ci]]),score=gene.data[[ci]]$t)
      cmap.drug.rank = drug.ref.profiles$drug.rank.matrix
      e1 = simpleError("Did not find the column named 'geneSymbol' in query data that contains the gene symbols in it.")
      e2 = simpleError("Did not find the column named 'score' in query data that contains the test statistics or any values that you would like to rank the genes.")

      cat("\n")
      cat("\n")
      message("Data preprocessing ...\n")
      cat("\n")
      if("score" %in% colnames(query.data)){
        if("geneSymbol" %in% colnames(query.data)){
          tmp = data_preprocess(query.data, cmap.drug.rank,connectivity = connectivity)
          query.data = tmp[[1]]
          cmap.drug.rank = tmp[[2]]
          rm(tmp)
        } else{
          stop(e1)
        }
      } else{
        stop(e2)
      }

      message("Identifying drug instance CEGs...\n")
      cat("\n")
      p_min = get_gene_pval('min',cmap.drug.rank,query.data)
      p_max = get_gene_pval('max',cmap.drug.rank,query.data)

      ##Select the smallest p value (between 2 p values) as the p value of the gene
      p_score = pmin(p_min,p_max)
      z_score = qnorm(p_score,lower.tail = F)
      CEG.pvals = get_CEGs(p_min, p_max, z_score,threshold = CEG.threshold)

      message("Calculating drug connectivity p values ...\n")
      cat("\n")
      drug.info = drug.ref.profiles$drug.info
      if(repurposing.unit == "drug"){
        drug.info$drug = drug.info$cmap_name
      } else if(repurposing.unit == "treatment"){
        drug.info$drug = drug.info$treatment
      } else{
        stop(simpleError("Please set the repurposing unit to either 'drug' or 'treatment'."))
      }

      drug.pvals = get_drug_pval(CEGsum = CEG.pvals$CEG.sumz.scores,drug.info = drug.info)

      drug.pvals = drug.pvals[order(drug.pvals$pval),]
      drugs = rownames(drug.pvals)
      drug.pvals$Drug.name = gsub("_.*","",drugs)
      drug.pvals$Drug.id = gsub(".*_","",drugs)
      rownames(drug.pvals) = NULL
      drug.pvals = drug.pvals[,c(3,4,1)]
      drug.pvals$FDR = p.adjust(drug.pvals$pval,method = "fdr")
      colnames(drug.pvals)[3] = "P.value"

      res = list(drug.pvals,drug.info,CEG.pvals$CEG.pvals)
      names(res) = c("drug.pvals","drug.info","CEG.pvals")
      res.list[[ci]] <- drug.pvals
  }
  names(res.list) <- names(gene.data)
  return(res.list)
}
