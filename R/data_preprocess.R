
##Match common genes between users query data and cmap
data_preprocess = function(query.data,cmap.drug.rank, connectivity){
  common.genes = intersect(query.data$geneSymbol,rownames(cmap.drug.rank))
  rownames(query.data) = query.data$geneSymbol
  query.data = query.data[common.genes,]
  if(connectivity == "negative"){
    ##Rank query data gene statistic scores from smallest to largest, opposite to cmap gene rank
    query.data$geneRank = rank(query.data$score,ties.method = "first")
  } else if(connectivity == "positive"){
    ##Rank query data gene statistic scores from largest to smallest, same with cmap gene rank
    query.data$geneRank = rank(-(query.data$score),ties.method = "first")
  }

  cmap.drug.rank = cmap.drug.rank[common.genes,]
  ##Re-rank drug rank matrix after excluding uncommon genes
  for(i in 1:ncol(cmap.drug.rank)){
    cmap.drug.rank[,i] = rank(cmap.drug.rank[,i])
  }
  return(list(query.data,cmap.drug.rank))
}

