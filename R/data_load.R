#' @export 
load_tnbc <- function(cache_dir=NULL) {
    # Download the data
    if (is.null(cache_dir)) {
        cache_dir <- BiocFileCache::getBFCOption("CACHE")
    }
    bfc <- BiocFileCache::BiocFileCache(cache_dir, ask=FALSE)

    url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE113197&format=file"
    gse113197_path <- BiocFileCache::bfcrpath(bfc, "gse113197", fpath=url)
    gse113197_temp <- tempfile()
    untar(gse113197_path, exdir=gse113197_temp)

    url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE123926&format=file"
    gse123926_path <- BiocFileCache::bfcrpath(bfc, "gse123926", fpath=url)
    gse123926_temp <- tempfile()
    untar(gse123926_path, exdir=gse123926_temp)

    # Load control samples
    sample_files <- list(
        Ind5 = "GSM3099847_Ind5_Expression_Matrix.txt.gz", 
        Ind6 = "GSM3099848_Ind6_Expression_Matrix.txt.gz",
        Ind7 = "GSM3099849_Ind7_Expression_Matrix.txt.gz"
    )
    sobjects <- list()

    for (i in seq_along(sample_files)) {
        sample_file <- gzfile(file.path(gse113197_temp, sample_files[[i]]))
        sample_name <- names(sample_files)[i]

        # Read the data and subset to the epithelial cells
        counts <- read.table(file = sample_file, header = TRUE, check.names = FALSE)
        row.names(counts) <- counts[, 1]
        counts <- counts[, -1]
        common <- intersect(colnames(counts), tnbc_cells[[sample_name]])
        counts <- counts[, common]

        meta <- data.frame(row.names = colnames(counts))
        meta$sample <- sample_name
        meta$type <- "control"
        sobjects[[sample_name]] <- SeuratObject::CreateSeuratObject(
            counts = counts, min.cells = 3, min.features = 200, meta.data = meta
        )
    }

    # Load case samples
    sample_files <- list(PDX110 = "GSM3516947_PDX110", PDX322 = "GSM3516948_PDX322")

    for (i in seq_along(sample_files)) {
        sample_name <- names(sample_files)[i]

        mtx_file = paste(sample_files[[i]], "matrix.mtx.gz", sep="-")
        cells_file = paste(sample_files[[i]], "barcodes.tsv.gz", sep="-")
        features_file = paste(sample_files[[i]], "genes.tsv.gz", sep="-")
        counts <- Seurat::ReadMtx(
            mtx = file.path(gse123926_temp, mtx_file),
            cells = file.path(gse123926_temp, cells_file), 
            features = file.path(gse123926_temp, features_file) 
        )

        meta <- data.frame(row.names = colnames(counts))
        meta$sample <- sample_name
        meta$type <- "control"
        sobjects[[sample_name]] <- SeuratObject::CreateSeuratObject(
            counts = counts, min.cells = 3, min.features = 200, meta.data = meta
        )
    }

    return(sobjects)
}