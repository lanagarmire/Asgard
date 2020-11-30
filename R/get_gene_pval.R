
##Beta test for gene p values
get_gene_pval = function(order, cmap.drug.rank,query.data){
  geneRank = as.matrix(replicate(ncol(cmap.drug.rank),query.data$geneRank))
  ##Keep min(x,y) to find the bottom ranked genes
  if (order == 'min'){
    order_stat = pmin(as.matrix(cmap.drug.rank),geneRank)
    p_val = 1 - pbeta((order_stat-1)/nrow(order_stat),1,2,lower.tail = T)
  }
  else if(order == 'max'){
    order_stat = pmax(as.matrix(cmap.drug.rank),geneRank)
    p_val = pbeta(order_stat/nrow(order_stat),2,1,lower.tail = T)
  }
  return(p_val)
}

