#' @title Drug Score
#' @description  The drug score is a comprehensive estimation of drug therapeutic effects acrossing all or selected single cell clusters. 
#' @details This function evaluates treatment efficacy and ranks drugs using therapeutics score, which integrates gene responses to multiple drugs, the proportion of genes, and cells treated by drugs.
#' @param SC.integrated A Seurat object of aligned single cells from Seurat.
#' @param Gene.data A list of differnential gene expression profiles for every cell type. It's from GetGene function.
#' @param Cell.type Select which clusters (cell types) to be used for drug score estimation. By default, it uses all clusters.
#' @param Drug.data A list of mono-drugs for every cluster. It's from GetDrug function.
#' @param FDA.drug.only logical; if TRUE, will only return FDA-approved drugs, else, will return all inputted drugs/compounds
#' @param GSE92742.gctx The gctx file contains drug responses from GSE92742 dataset (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92742).
#' @param GSE70138.gctx The gctx file contains drug responses from GSE70138 dataset (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70138).
#' @param Case A vector contains case sample names.
#' @param Tissue Reference tissue. If one used lung_rankMatrix.txt in GetDrugRef function, then the Reference tissue is lung.
#' @return A data frame of drug score, P-value and FDR.
#' @export
#' @import cmapR


DrugScore <- function(SC.integrated=SC.data,
                            Gene.data=Gene.list,
                            Cell.type=NULL,
                            Drug.data=Drug.ident.res,
                            FDA.drug.only=TRUE,
                            GSE92742.gctx=NULL,
                            GSE70138.gctx=NULL,
                            Case=NULL,
                            Tissue="breast"
){
    if(length(Cell.type)>0){
      Cell.type=intersect(Cell.type, unique(SC.integrated@meta.data$celltype))
      SC.integrated=subset(SC.integrated, celltype %in% Cell.type)
      Drug.data=Drug.data[Cell.type]
      Gene.data=Gene.data[Cell.type]
    }
    ##Cell proportion
    cells <- SC.integrated@meta.data
    if(length(Case)>0){
      cells <- subset(cells,sample %in% Case)
      SC.integrated=subset(SC.integrated, sample %in% Case)
    }
    cells <- cells$celltype
    cell.count <- table(cells)
    cell.count <- cell.count[which(cell.count>3)]
    cells.freq <- round(100*cell.count/length(cells),2)
    cells.freq.rank <- rank(as.vector(cells.freq))

    ##Load drug data
    Drug.list <- data.frame()
    for(i in names(Drug.data)){
      Cd <- Drug.data[[i]]
      Cd <- Cd[!duplicated(Cd$Drug.name),]
      Drugs <- Cd$Drug.name
      if(FDA.drug.only==TRUE){
      Drugs <- intersect(Drugs,FDA.drug)
      }
      if(length(Drugs)>0){
      Cd <- subset(Cd, Drug.name %in% Drugs)
      FDRs <- Cd$FDR
      Pvalue <- Cd$P.value
      temp <- data.frame(Drug=Drugs,Cluster=i,Size=cells.freq[i],P.value=Pvalue,FDR=FDRs,row.names = NULL)
      Drug.list <- rbind(Drug.list,temp)
      }
    }
    Drug.list <- unique(Drug.list)
    Drug.list$P.value <- Drug.list$P.value
    Drug.list$w.size <- Drug.list$Size*(-log10(Drug.list$FDR))
    Drug.list[is.na(Drug.list)] <- 0
    Drug.coverage <- tapply(Drug.list$w.size, Drug.list$Drug,sum)
    C.Drugs <- rownames(Drug.coverage)
    
    ##Combine P value
    if(length(unique(names(Drug.data)))>1){
       Combined.Pvalue <- tapply(Drug.list$P.value, Drug.list$Drug, CombineP)
    }else{
      Combined.Pvalue <- Drug.list$P.value
      names(Combined.Pvalue) <- Drug.list$Drug
    }
  
    ##Cell line information
    cells <- subset(cell_data,primary_site == Tissue)$cell_id

    ##Load experiment information
    data_infor1 <- col_meta_GSE92742[,c("sig_id","pert_iname")]
    row.names(data_infor1) <- data_infor1$sig_id
    idx <- which(col_meta_GSE92742$cell_id %in% cells & col_meta_GSE92742$pert_iname %in% C.Drugs)
    sig_ids <- col_meta_GSE92742$sig_id[idx]
    data_infor1 <- data_infor1[sig_ids,]

    ##Load drug response
    my_ds <- parse_gctx(GSE92742.gctx, cid=sig_ids)
    gene.data <- as.data.frame(my_ds@mat)
    gene.data$geneid <- row.names(gene.data)
    treatments <- colnames(gene.data)
    treatments <- setdiff(treatments,"geneid")
    data <- merge(gene.data,gene_meta,by.x="geneid",by.y="pr_gene_id")
    data1 <- data[,c("pr_gene_symbol",treatments)]

    ##Load experiment information
    data_infor2 <- col_meta_GSE70138[,c("sig_id","pert_iname")]
    row.names(data_infor2) <- data_infor2$sig_id
    idx <- which(col_meta_GSE70138$cell_id %in% cells & col_meta_GSE70138$pert_iname %in% C.Drugs)
    sig_ids <- col_meta_GSE70138$sig_id[idx]
    data_infor2 <- data_infor2[sig_ids,]

    ##Load drug response
    sig_ids <- col_meta_GSE70138$sig_id[idx]
    my_ds <- parse_gctx(GSE70138.gctx, cid=sig_ids)
    gene.data <- as.data.frame(my_ds@mat)
    gene.data$geneid <- row.names(gene.data)
    treatments <- colnames(gene.data)
    treatments <- setdiff(treatments,"geneid")
    data <- merge(gene.data,gene_meta,by.x="geneid",by.y="pr_gene_id")
    data2 <- data[,c("pr_gene_symbol",treatments)]
    data <- merge(data1,data2,by="pr_gene_symbol")
    row.names(data) <- data[,1]
    data <- data[,-1]
    data_infor <- rbind(data_infor1,data_infor2)

    ##Drug score
    D.genes <- list()
    for(i in names(Gene.data)){
      Cd <- Gene.data[[i]]
      Cd <- subset(Cd, adj.P.Val<0.05)
      D.genes.temp <- list(temp=rownames(Cd))
      D.genes <- cbind(D.genes,D.genes.temp)
    }
    D.genes <- Reduce(intersect,D.genes)
    Gene.expression <- data.frame()
    for(i in names(Gene.data)){
      Cd <- Gene.data[[i]]
      if(nrow(Gene.expression)==0){
        Gene.expression <- data.frame(Score=Cd[D.genes,"score"])
      }else{
      Gene.expression.temp <- data.frame(Score=Cd[D.genes,"score"])
      Gene.expression <- cbind(Gene.expression,Gene.expression.temp)
      }
    }
    Gene.expression <- as.data.frame(Gene.expression)
    Gene.expression <- as.matrix(Gene.expression)
    row.names(Gene.expression) <- D.genes
    D.gene.expression <- apply(Gene.expression,1,mean)
    names(D.gene.expression) <- D.genes
    Single.treated.score.list <- NULL
    for(Drug in C.Drugs){
      D.genes.treated <- NULL
      drug.treatments <- subset(data_infor,pert_iname == Drug)$sig_id
      if(length(drug.treatments)>1){
        drug.responses <- data[,drug.treatments]
        drug.responses.mean <- apply(drug.responses,1,mean)
      }else{
        drug.responses <- data[,drug.treatments]
        drug.responses.mean <- drug.responses
      }
      D.D.genes <- intersect(names(D.gene.expression),names(drug.responses.mean))
      D.genes.treated <- -D.gene.expression[D.D.genes]*drug.responses.mean[D.D.genes]
      D.genes.treated <- D.genes.treated[which(D.genes.treated>0)]
      Mean.treated <- mean(D.genes.treated)
      Ratio.treated <- length(D.genes.treated)/length(D.D.genes)
      Coverage.treated <- Drug.coverage[Drug]/100
      Treated.score <- (Ratio.treated*Coverage.treated)
      Single.treated.score.list <- c(Single.treated.score.list,Treated.score)
    }
    Res.table <- data.frame(Drug.therapeutic.score=Single.treated.score.list[C.Drugs],P.value=Combined.Pvalue[C.Drugs],FDR=p.adjust(Combined.Pvalue[C.Drugs],method = "BH"))
    return(Res.table)

}
