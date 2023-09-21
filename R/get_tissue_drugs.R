#' @export 
get_tissue_drugs <- function(
    tissue, cell_info, gene_info, gse70138_sig_info, gse92742_sig_info, 
    gse70138_gctx, gse92742_gctx, output_dir
) {    
    # Get cell line information for L1000
    cell_meta <- read.table(cell_info, sep="\t", header=TRUE, quote="")
    available_tissues <- unique(cell_meta$primary_site)
    available_tissues <- available_tissues[which(available_tissues != "-666")]

    # Check if the tissue that needs to be processed is available
    stopifnot(tissue %in% available_tissues)

    # Get cell lines names
    cell_ids <- which(cell_meta$primary_site == tissue)
	cell_names <- cell_meta$cell_id[cell_ids]

    drugs1 <- get_series_drugs(tissue, cell_names, gse70138_sig_info, 
                               gse70138_gctx)
    drugs2 <- get_series_drugs(tissue, cell_names, gse92742_sig_info, 
                               gse92742_gctx)

    # Combine drugs from two GSE series
    expressions = list()
    sig_meta = list()
    i <- 1
    if (!is.null(drugs1)) {
        expressions[[i]] <- drugs1$expressions
        sig_meta[[i]] <- drugs1$sig_meta
        i <- i+1
    }

    if (!is.null(drugs2)) {
        expressions[[i]] <- drugs2$expressions
        sig_meta[[i]] <- drugs2$sig_meta
    }

    if (length(expressions) == 0) {
        stop(sprintf("Warning: no drugs found for %s.", tissue))
    } else {
        expressions <- dplyr::bind_cols(expressions)
        sig_meta <- dplyr::bind_rows(sig_meta)
    }

    # Get gene metadata
    gene_meta <- read.table(gene_info, sep="\t", header=TRUE, quote = "")

    expressions <- expressions[gene_meta$pr_gene_id, ]
    rownames(expressions) <- gene_meta$pr_gene_symbol

    return(list(
        expressions=expressions, gene_meta=gene_meta, sig_meta=sig_meta
    ))
}

get_series_drugs <- function(tissue, cell_names, sig_info, gctx) {
    # Get signature information
    sig_meta <- read.table(sig_info, sep="\t", header=TRUE, quote="")

    # Get signatures related to the tissue we are analyzing
    if(tissue == "breast") {
        idx <- which(sig_meta$cell_id %in% cell_names & 
                     sig_meta$pert_type == "trt_cp" & 
                     sig_meta$pert_id != "BRD-K18910433") #TODO ??
    } else {
        idx <- which(sig_meta$cell_id %in% cell_names & 
                     sig_meta$pert_type == "trt_cp")
    }
    sig_ids <- sig_meta$sig_id[idx]

    # Remove repeating signature #TODO: ??
    rm_ids <- grep('REP\\.', sig_ids)
    if(length(rm_ids)>0){
        sig_ids <- sig_ids[-rm_ids]
    }

    # Get gene-signature interactions
    if (length(sig_ids) > 0) {
        exprs <- cmapR::parse_gctx(gctx, cid = sig_ids)@mat
        
        sig_meta <- subset(sig_meta, sig_id %in% sig_ids)
        sig_meta <- restruct_sig_meta(sig_meta)

        out = list(expressions = as.data.frame(exprs), sig_meta = sig_meta)
        return(out)
    } else {
        return(invisible(NULL))
    }
}

restruct_sig_meta <- function(df) {
    return(data.frame(
        sig_id = df$sig_id,
        cmap_name = paste(df$pert_iname, df$pert_id, sep="_"),
        concentration = df$pert_idose,
        duration = df$pert_itime,
        cell_id = df$cell_id,
        catalog_name = df$pert_id,
        treatment = paste(df$pert_iname, df$sig_id, sep="_")
    ))
}