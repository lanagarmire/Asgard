
##Get drug p value based on CEGscore
get_drug_pval = function(CEGsum,drug.info){
  treat_drug_ks = matrix(0,ncol=1,nrow = length(unique(drug.info$drug)))
  treat_drug_ks = as.data.frame(treat_drug_ks)
  rownames(treat_drug_ks) = unique(drug.info$drug)
  colnames(treat_drug_ks) = 'pval'

  for(i in 1:nrow(treat_drug_ks)){
    indiv_drug = drug.info[(drug.info$drug == rownames(treat_drug_ks)[i]),]$instance_id
    indiv_drug  <-  as.character(indiv_drug)
    indiv_drug_score = CEGsum[indiv_drug]
    rest_score = CEGsum[setdiff(names(CEGsum),indiv_drug)]

    ##k-s test: one drug drug.info v.s. other drug.info
    options(warn = -1)
    treat_drug_ks$pval[i] = (ks.test(indiv_drug_score,rest_score,alternative = 'less'))$p.value
  }
  treat_drug_ks$drug = sapply(rownames(treat_drug_ks),function(x){strsplit(x,split="_")[[1]][1]})

  return(treat_drug_ks)
}

