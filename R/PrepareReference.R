#' @title Prepare Drug Reference.
#' @description  Prepare tissue specific drug reference Profiles from L1000 drug response data.
#' @details This function converts L1000 data to the tissue specific drug rank matrix.
#' @param cell.info The local path and the name of the cell.info text file. It's downloaded from https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/GSE70138_Broad_LINCS_cell_info_2017-04-28.txt.gz .
#' @param gene.info The local path and the name of the gene.info text file. It's downloaded from https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/GSE70138_Broad_LINCS_gene_info_2017-03-06.txt.gz .
#' @param GSE70138.sig.info The local path and the name of the cell.info text file. It's downloaded from https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/GSE70138_Broad_LINCS_sig_info_2017-03-06.txt.gz .
#' @param GSE92742.sig.info The local path and the name of the cell.info text file. It's downloaded from https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_sig_info.txt.gz .
#' @param GSE70138.gctx The local path and the name of the cell.info text file. It's downloaded from https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx.gz .
#' @param GSE92742.gctx The local path and the name of the cell.info text file. It's downloaded from https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx.gz .
#' @param Output.Dir The output directory for the generated files.
#' @export
#' @import cmapR

PrepareReference <- function(cell.info = NULL,
                    gene.info = NULL,
                    GSE70138.sig.info = NULL,
                    GSE92742.sig.info = NULL,
                    GSE70138.gctx = NULL,
                    GSE92742.gctx = NULL,
                    Output.Dir = "./"){
      cell_data<-read.table(file=cell.info,sep="\t",header = T,quote = "")
      tissues<-unique(as.character(cell_data$primary_site))
      tissues<-tissues[which(tissues!="-666")]
      for (tissue in tissues){
          print(tissue)
          cell_data<-read.table(file=cell.info,sep="\t",header = T,quote = "")
          cell_ids<-which(cell_data$primary_site == tissue)
          cell_names <- cell_data$cell_id[cell_ids]
          ds_path <- GSE70138.gctx
          col_meta_path <- GSE70138.sig.info
          col_meta <- read.delim(col_meta_path, sep="\t", stringsAsFactors=F)
          if(tissue == "breast"){
          idx <- which(col_meta$cell_id %in% cell_names & col_meta$pert_type == "trt_cp" & col_meta$pert_id!="BRD-K18910433")
          }else{
          idx <- which(col_meta$cell_id %in% cell_names & col_meta$pert_type == "trt_cp")
          }
          sig_ids <- col_meta$sig_id[idx]
          length1<-length(sig_ids)
          if(length1 > 0){
          my_ds <- parse_gctx(ds_path, cid=sig_ids)
          myrank <- function(x){
            temp<-rank(-x,ties.method ="min")
            return(temp)
          }
          rank_matrix1<-apply(my_ds@mat,2,myrank)
          rank_matrix1<-as.data.frame(rank_matrix1)
          }
          cell_data<-read.table(file=cell.info,sep="\t",header = T,quote = "")
          cell_ids<-which(cell_data$primary_site == tissue)
          cell_names <- cell_data$cell_id[cell_ids]
          ds_path <- GSE92742.gctx
          col_meta_path <- GSE92742.sig.info
          col_meta <- read.delim(col_meta_path, sep="\t", stringsAsFactors=F)
          if(tissue == "breast"){
          idx <- which(col_meta$cell_id %in% cell_names & col_meta$pert_type == "trt_cp" & col_meta$pert_id!="BRD-K18910433")
          }else{
          idx <- which(col_meta$cell_id %in% cell_names & col_meta$pert_type == "trt_cp")
          }
          sig_ids <- col_meta$sig_id[idx]
          length2<-length(sig_ids)
          if(length2 > 0){
          my_ds <- parse_gctx(ds_path, cid=sig_ids)
          myrank <- function(x){
            temp<-rank(-x,ties.method ="min")
            return(temp)
           }
          rank_matrix2<-apply(my_ds@mat,2,myrank)
          rank_matrix2<-as.data.frame(rank_matrix2)
          }

          if(length1 > 0 & length2 > 0){
          rank_matrix<-cbind(rank_matrix1,rank_matrix2)
          }else if(length1 > 0 & length2 == 0){
          rank_matrix<-rank_matrix1
          }else if(length1 == 0 & length2 > 0){
          rank_matrix<-rank_matrix2
          }
          if(length1 > 0 |  length2 > 0){
              colnames(rank_matrix)<-gsub(":","_",colnames(rank_matrix))
              cnames<-colnames(rank_matrix)
              colnames(rank_matrix)<-1:length(cnames)
              dcnames<-colnames(rank_matrix)
              rank_matrix$probe_id<-row.names(rank_matrix)
              rank_matrix <- rank_matrix[,c('probe_id', dcnames)]
              filename<-paste(Output.Dir,gsub(" ","-",tissue),"_rankMatrix.txt",sep = "")
              write.table(rank_matrix,file=filename,quote=FALSE,row.names = FALSE,sep = "\t")

              gene_data<-read.table(file=gene.info,sep="\t",header = T,quote = "")
              my_gene_info<-gene_data[,1:2]
              colnames(my_gene_info)<-c("ID","Gene.Symbol")
              filename<-paste(Output.Dir,gsub(" ","-",tissue),"_gene_info.txt",sep = "")
              write.table(my_gene_info,file=filename,quote=FALSE,row.names = FALSE,sep = "\t")

              sig_data<-read.table(file=GSE70138.sig.info,sep="\t",header = T,quote = "")
              sig_data$sig_id<-gsub(":","_",sig_data$sig_id)
              my_drug_info<-data.frame(instance_id=sig_data$sig_id,cmap_name=paste(sig_data$pert_iname,sig_data$pert_id,sep="_"),concentration..M=sig_data$pert_idose,duration..h=sig_data$pert_itime,cell2=sig_data$cell_id,catalog_name=sig_data$pert_id,treatment=paste(sig_data$pert_iname,"_",sig_data$sig_id,sep = ""))
              my_drug_info1<-subset(my_drug_info,instance_id %in% cnames)

              sig_data<-read.table(file=GSE92742.sig.info,sep="\t",header = T,quote = "")
              sig_data$sig_id<-gsub(":","_",sig_data$sig_id)
              my_drug_info<-data.frame(instance_id=sig_data$sig_id,cmap_name=paste(sig_data$pert_iname,sig_data$pert_id,sep="_"),concentration..M=sig_data$pert_idose,duration..h=sig_data$pert_itime,cell2=sig_data$cell_id,catalog_name=sig_data$pert_id,treatment=paste(sig_data$pert_iname,"_",sig_data$sig_id,sep = ""))
              my_drug_info2<-subset(my_drug_info,instance_id %in% cnames)

              my_drug_info<-rbind(my_drug_info1,my_drug_info2)
              my_drug_info$instance_id<-1:length(my_drug_info$instance_id)
              filename<-paste(Output.Dir,gsub(" ","-",tissue),"_drug_info.txt",sep = "")
              write.table(my_drug_info,file=filename,quote=FALSE,row.names = FALSE,sep = "\t")
          }

    }
}
