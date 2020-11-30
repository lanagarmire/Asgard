
##Get CEGs and CEG's sumz scores
get_CEGs = function(p_min, p_max, z_score,threshold){
  CEG.pvals = list()
  CEG.pvals$down = p_min
  CEG.pvals$up = p_max
  CEGz= numeric(ncol(z_score))
  for(i in 1:ncol(z_score)){
    CEGz[i] = sum(z_score[which(z_score[,i] >= qnorm(threshold,lower.tail = F)),i])
  }
  names(CEGz) = colnames(z_score)
  res = list(CEGz, CEG.pvals)
  names(res) = c("CEG.sumz.scores","CEG.pvals")
  return(res)
}

