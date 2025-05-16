#' @title Calculate drug score
#' @description The drug score is a comprehensive estimation of drug therapeutic 
#' effects using all or a selected set of clusters. 
#' @details This function calculates drug score using cellular proportion of 
#' clusters, the significance of reversal in DEGs' expressions, and the ratio of 
#' the reversed genes. 
#' @param cell_metadata A data.frame of cell metadata. It must have a column 
#' named 'cluster' indicating which cluster cells belong, and a column named 
#' 'sample' indicating which sample cells belong. 
#' @param cluster_degs A list of differential gene expression profiles for 
#' each cluster.
#' @param cluster_drugs Drug repurposing result from GetDrug function.
#' @param tissue Reference tissue. If one used 'lung_rankMatrix.txt' in 
#' GetDrugRef function, then the Reference tissue is lung. Please use " " 
#' instead of "-" in tissue name. For example, while 
#' 'haematopoietic-and-lymphoid-tissue' is the prefix of the drug reference 
#' files, the corresponding tissue name is "haematopoietic and lymphoid tissue".
#' @param gse70138_gctx_path The gctx file contains drug responses from GSE70138 
#' dataset (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70138).
#' @param gse92742_gctx_path The gctx file contains drug responses from GSE92742 
#' dataset (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92742)..
#' @param clusters Select which clusters (cell types) to be used for drug score 
#' estimation. By default, it uses all clusters.
#' @param case A vector containing case sample names.
#' @param fda_drugs_only logical; if TRUE, will only return FDA-approved drugs, 
#' else, will return all drugs/compounds.
#' @return A data frame of drug score, P-value and FDR.
#' @export
#' @import cmapR
DrugScore <- function(cell_metadata, cluster_degs, cluster_drugs, tissue,
					  gse70138_gctx_path, gse92742_gctx_path, 
					  clusters = NULL, case = NULL, fda_drugs_only = TRUE) {

	# Subset input data to the set of clusters we are interested in 
    if (length(clusters) > 0) {
    	clusters = intersect(clusters, unique(cell_metadata$cluster))
      	cell_metadata = subset(cell_metadata, cluster %in% clusters)
      	cluster_drugs = cluster_drugs[clusters]
      	cluster_degs = cluster_degs[clusters]
    }

    # Calculate cluster proportions in diseased tissue
    if (length(case) > 0) {
      	cell_metadata <- subset(cell_metadata, sample %in% case)
    }
    clustering <- cell_metadata$cluster
    cluster_sizes <- table(clustering)
    cluster_sizes <- cluster_sizes[which(cluster_sizes > 3)]
    cluster_prop <- round(100*cluster_sizes/nrow(cell_metadata), 2) 

    # Combine cluster drugs into a single data frame
    drug_list <- data.frame()
    for (i in names(cluster_drugs)) {
    	ith_cluster_drugs <- cluster_drugs[[i]]
		drug_names <- ith_cluster_drugs$Drug.name
      	ith_cluster_drugs <- ith_cluster_drugs[!duplicated(drug_names), ]

		# Subset to FDA drugs
      	if (fda_drugs_only) {
      		drug_names <- intersect(drug_names, FDA.drug)
      	}

      	if (length(drug_names)>0) {
      		ith_cluster_drugs <- subset(ith_cluster_drugs, Drug.name %in% drug_names)
      		fdrs <- ith_cluster_drugs$FDR
      		p_values <- ith_cluster_drugs$P.value
      		
			temp <- data.frame(
				drug = drug_names, 
				cluster = i,
				cluster_prop = cluster_prop[i],
				p_value = p_values,
				fdr = fdrs,
				row.names = NULL
			)
      		drug_list <- rbind(drug_list, temp)
      	}
    }
    drug_list <- unique(drug_list)
    drug_list$weighted_prop <- drug_list$cluster_prop*(-log10(drug_list$fdr))
    drug_list[is.na(drug_list)] <- 0

    drug_coverage <- tapply(drug_list$weighted_prop, drug_list$drug, sum)
    drugs <- rownames(drug_coverage)

    # Combine cluster spesific p-values of drugs
    if(length(unique(names(cluster_drugs)))>1){
       	combined_p_values <- tapply(drug_list$p_value, drug_list$drug, CombineP)
    }else{
      	combined_p_values <- drug_list$p_value
      	names(combined_p_values) <- drug_list$drug
    }
  
	# Cell line information
    cell_lines <- subset(cell_data, primary_site == tissue)$cell_id

    # Load drugs metadata for GSE92742 and subset it to tissue of interest and 
	# drugs of interest
    drug_metadata_92742 <- col_meta_GSE92742[, c("sig_id", "pert_iname")]
    row.names(drug_metadata_92742) <- drug_metadata_92742$sig_id
    idx <- which(col_meta_GSE92742$cell_id %in% cell_lines & 
				 col_meta_GSE92742$pert_iname %in% drugs)
    sig_ids <- col_meta_GSE92742$sig_id[idx]
    drug_metadata_92742 <- drug_metadata_92742[sig_ids, ]

    # Load drug response for GSE92742
    exprs <- as.data.frame(parse_gctx(gse92742_gctx_path, cid=sig_ids)@mat)
    treatments <- colnames(exprs)
	exprs$gene_id <- row.names(exprs)
    tmp <- merge(exprs, gene_meta, by.x="gene_id", by.y="pr_gene_id")
    drug_responses_92742 <- tmp[, c("pr_gene_symbol", treatments)]

    # Load drugs metadata for GSE70138 and subset it to tissue of interest and 
	# drugs of interest
    drug_metadata_70138 <- col_meta_GSE70138[, c("sig_id", "pert_iname")]
    row.names(drug_metadata_70138) <- drug_metadata_70138$sig_id
    idx <- which(col_meta_GSE70138$cell_id %in% cell_lines & 
				 col_meta_GSE70138$pert_iname %in% drugs)
    sig_ids <- col_meta_GSE70138$sig_id[idx]
    drug_metadata_70138 <- drug_metadata_70138[sig_ids, ]

    # Load drug response for GSE70138
    exprs <- as.data.frame(parse_gctx(gse70138_gctx_path, cid=sig_ids)@mat)
    treatments <- colnames(exprs)
    exprs$gene_id <- row.names(exprs)
	tmp <- merge(exprs, gene_meta, by.x="gene_id", by.y="pr_gene_id")
    drug_responses_70138 <- tmp[, c("pr_gene_symbol", treatments)]

    drug_responses <- merge(drug_responses_92742, drug_responses_70138, 
							by="pr_gene_symbol")
    row.names(drug_responses) <- drug_responses[, 1]
    drug_responses <- drug_responses[, -1]
    drug_metadata <- rbind(drug_metadata_92742, drug_metadata_70138)

	# Find DEGs that are common to all clusters
    common_degs <- list()
    for (i in names(cluster_degs)) {
    	ith_cluster_degs <- cluster_degs[[i]]
      	ith_cluster_degs <- subset(ith_cluster_degs, adj.P.Val < 0.05)
		if (length(ith_cluster_degs) > 0) {
	    	common_degs[[i]] <- rownames(ith_cluster_degs)
		}
    }
	common_degs <- Reduce(intersect, common_degs)

	# Combine cluster specific DEG scores into a matrix
	deg_scores <- data.frame()
    for (i in names(cluster_degs)) {
    	ith_cluster_degs <- cluster_degs[[i]]
    	if (nrow(deg_scores) == 0) {
    		deg_scores <- data.frame(score = ith_cluster_degs[common_degs, "score"])
    	} else {
    	    tmp <- data.frame(score = ith_cluster_degs[common_degs,"score"])
        	deg_scores <- cbind(deg_scores, tmp)
       }
    }
    deg_scores <- as.matrix(deg_scores)
    row.names(deg_scores) <- common_degs

    deg_scores_mean <- apply(deg_scores, 1, mean)
    names(deg_scores_mean) <- common_degs

	# Calculate drug score
    drug_scores <- list()
    for (drug in drugs) {
		# Get response from CMap
		treatments <- subset(drug_metadata, pert_iname == drug)$sig_id
		if (length(treatments) > 1) {
			curr_drug_response <- drug_responses[, treatments]
			mean_response <- apply(curr_drug_response, 1, mean)
		} else {
			curr_drug_response <- drug_responses[, treatments]
			mean_response <- curr_drug_response
		}

		drug_stats <- drug_list[drug_list$drug == drug, ]
		drug_score <- 0
		for (i in names(cluster_degs)) {
			cluster_prop <- drug_stats[drug_stats$cluster == i, "cluster_prop"]
			fdr <- drug_stats[drug_stats$cluster == i, "fdr"]
			p_value <- drug_stats[drug_stats$cluster == i, "p_value"]

			ith_cluster_degs <- cluster_degs[[i]]
      		ith_cluster_degs <- subset(ith_cluster_degs, adj.P.Val < 0.05)

			treatable_degs <- intersect(row.names(ith_cluster_degs), names(mean_response))
			if (length(treatable_degs > 0)) {
				deg_scores <- ith_cluster_degs[treatable_degs, "score"]

				treated_degs <- -deg_scores*mean_response[treatable_degs]
				treated_degs <- treated_degs[which(treated_degs > 0)]

				treated_degs_ratio <- length(treated_degs)/length(treatable_degs)
				drug_score <- drug_score +
					(cluster_prop/100)*(-log10(fdr))*treated_degs_ratio
			}
	    }
		
		drug_scores[[drug]] <- drug_score
    }
	drug_scores <- t(as.data.frame(drug_scores))

    out <- data.frame(
		Drug.therapeutic.score = drug_scores,
		P.value = combined_p_values[drugs],
		FDR = p.adjust(combined_p_values[drugs], method = "BH")
	)
    return(out)

}
