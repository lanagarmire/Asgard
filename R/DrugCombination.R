#' @title Treatment Efficacy of the Drug Combination.
#' @description  It evaluates treatment efficacy to identify drug combinations that can best reverse the target genes’ expression in diseased cells in case samples.
#' @details This function evaluates treatment efficacy and ranks drug combinations using therapeutics score, which integrates gene responses to multiple drugs, the proportion of genes, and cells treated by combined drugs.
#' @param SC.integrated A Seurat object of aligned single cells from SCalignment function.
#' @param Gene.data A list of differnential gene expression profiles for every cell type. It's from GetGene function.
#' @param Drug.data A list of mono-drugs for every cell type. It's from GetDrug function.
#' @param Drug.FDR The FDR threshold to select drug. The default value is 0.1.
#' @param FDA.drug.only logical; if TRUE, will only return FDA-approved drugs.
#' @param Combined.drugs The number of drugs in a combination. The default value is 2.
#' @param GSE92742.gctx The gctx file contains drug responses from GSE92742 dataset (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92742).
#' @param GSE70138.gctx The gctx file contains drug responses from GSE70138 dataset (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70138).
#' @param Case A vector contains names of case samples.
#' @param Tissue Reference tissue. If one used lung_rankMatrix.txt in GetDrugRef function, then the Reference tissue is lung.
#' @return A data frame of drug combinations with therapeutics scores and FDR.
#' @export
#' @import cmapR


DrugCombination <- function(SC.integrated=SC.data,
                            Gene.data=Gene.list,
                            Drug.data=Drug.ident.res,
                            Drug.FDR=0.1,
                            FDA.drug.only=TRUE,
                            Combined.drugs=2,
                            GSE92742.gctx=NULL,
                            GSE70138.gctx=NULL,
                            Case=NULL,
                            Tissue="breast"
){
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
      Cd <- Cd[!duplicated(Cd$Drug.name),]
      #Cd <- subset(Cd, FDR<Drug.FDR)
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
    Drug.list$w.size <- Drug.list$Size*(-log10(Drug.list$FDR))
    Drug.list[is.na(Drug.list)] <- 0
    Drug.coverage <- tapply(Drug.list$w.size, Drug.list$Drug,sum)
    raw.raw.Drug.list <- Drug.list
    Drug.list <- subset(Drug.list, FDR<Drug.FDR)
    Drug.combinations <- combn(unique(Drug.list$Drug),Combined.drugs)
    Select.combnation <- function(x){
      temp.list <- subset(raw.raw.Drug.list,Drug %in% x)
      #temp.list <- unique(temp.list[,2:3])
      temp.size <- sum(temp.list$w.size)
      return(temp.size)
    }
    label<-apply(Drug.combinations,2,Select.combnation)
    Selected.Drug.combinations <- Drug.combinations[,which(label>0)]
    Selected.Drug.combinations.coverage <- label[which(label>0)]
    C.Drugs <- unique(as.vector(Selected.Drug.combinations))

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

    ##Combination score
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
      drug.responses <- data[,drug.treatments]
      drug.responses.mean <- apply(drug.responses,1,mean)
      D.D.genes <- intersect(names(D.gene.expression),names(drug.responses.mean))
      D.genes.treated <- -D.gene.expression[D.D.genes]*drug.responses.mean[D.D.genes]
      D.genes.treated <- D.genes.treated[which(D.genes.treated>0)]
      D.genes.treated <- D.genes.treated
      Mean.treated <- mean(D.genes.treated)
      Ratio.treated <- length(D.genes.treated)/length(D.D.genes)
      Coverage.treated <- Drug.coverage[Drug]/100
      Treated.score <- (Ratio.treated*Coverage.treated)
      Single.treated.score.list <- c(Single.treated.score.list,Treated.score)
    }
    Combination.treated.score <- function(Drugs){
      D.genes.treated<-NULL
      for(drug in Drugs){
        drug.treatments <- subset(data_infor,pert_iname == drug)$sig_id
        drug.responses <- data[,drug.treatments]
        drug.responses.mean <- apply(drug.responses,1,mean)
        D.D.genes <- intersect(names(D.gene.expression),names(drug.responses.mean))
        D.genes.treated.temp <- -D.gene.expression[D.D.genes]*drug.responses.mean[D.D.genes]
        D.genes.treated <- cbind(D.genes.treated,D.genes.treated.temp)
      }
      remove <- which(rowSums(D.genes.treated<0)==length(Drugs))
      D.genes.combination <- D.genes.treated[-remove,]
      scores <- apply(D.genes.combination,1,mean)
      temp.scores <- scores
      return(temp.scores)
    }
    Score.list <- apply(Selected.Drug.combinations, 2, Combination.treated.score)
    Combination.treated.ratio <- function(Drugs){
      D.genes.treated<-NULL
      for(drug in Drugs){
        drug.treatments <- subset(data_infor,pert_iname == drug)$sig_id
        drug.responses <- data[,drug.treatments]
        drug.responses.mean <- apply(drug.responses,1,mean)
        D.D.genes <- intersect(names(D.gene.expression),names(drug.responses.mean))
        D.genes.treated.temp <- -D.gene.expression[D.D.genes]*drug.responses.mean[D.D.genes]
        D.genes.treated <- cbind(D.genes.treated,D.genes.treated.temp)
      }
      remove <- which(rowSums(D.genes.treated<0)==length(Drugs))
      D.genes.combination <- D.genes.treated[-remove,]
      scores <- apply(D.genes.combination,1,mean)
      Ratio.treated <- length(which(scores>0))/length(D.D.genes)
      temp.scores <- Ratio.treated
      return(temp.scores)
    }
    Ratio.list <- apply(Selected.Drug.combinations, 2, Combination.treated.ratio)
    ref.score <- unlist(Score.list)
  	P.value <- function(Score) {
  	  if(length(Score)>1 && length(ref.score)>1){
  	  temp <- ks.test(Score, ref.score)
  	  p.value <- temp$p.value
  	  return(p.value)
  	  }else{
  		return(1)
  	  }
  	}
    pvalues <- unlist(suppressWarnings(lapply(Score.list, P.value)))
    combination.scores <- unlist(suppressWarnings(lapply(Ratio.list,mean)))
    Combination.table <- as.data.frame(t(Selected.Drug.combinations))
    for(d in 1:Combined.drugs){
     Combination.table <- cbind(Combination.table, Single.treated.score.list[Combination.table[,d]])
    }
    neg.combination.scores <- which(combination.scores<0)
    combination.scores[neg.combination.scores] <- -combination.scores[neg.combination.scores]
    Combination.table$Combination.therapeutic.score <- (Selected.Drug.combinations.coverage*combination.scores/100)
    Combination.table$Combination.therapeutic.score[neg.combination.scores] <- -Combination.table$Combination.therapeutic.score[neg.combination.scores]
    Combination.table$P.value <- pvalues
    Combination.table$FDR <- p.adjust(pvalues, method = "BH")
    colnames(Combination.table)[1:Combined.drugs] <- paste0("Drug",1:Combined.drugs)
    colnames(Combination.table)[(Combined.drugs+1):(2*Combined.drugs)] <- paste0("Drug",1:Combined.drugs,".therapeutic.score")
    return(Combination.table)
}
