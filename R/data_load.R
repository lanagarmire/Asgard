#' Load Triple Negative Breast Cancer dataset analyzed in the Asgard publication
#' 
#' This function load breast cancer datasets studied in the paper. At the first 
#' run, the function will downlaod necessary datasets from GSE and save them 
#' under a cache directory using `BiocFileCache`` library. 
#' 
#' @param cache_dir The directory to save downloaded datasets. If not giving 
#' the datasets are saved under the default cache directory of `BiocFileCache` 
#' library. 
#' 
#' @return A list of Seurat objects.  
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

    # Remove temp files 
    unlink(gse113197_temp, recursive=TRUE)
    unlink(gse123926_temp, recursive=TRUE)

    return(sobjects)
}

#' Load Luekemia dataset analyzed in the Asgard publication
#' 
#' This function load Luekimia datasets studied in the paper. At the first 
#' run, the function will downlaod necessary datasets from GSE and save them 
#' under a cache directory using `BiocFileCache`` library. 
#' 
#' @param cache_dir The directory to save downloaded datasets. If not giving,
#' the datasets are saved under the default cache directory of `BiocFileCache` 
#' library. 
#' 
#' @return A list of Seurat objects.  
#' @export 
load_luekemia <- function(cache_dir=NULL) {
    # Download the data
    if (is.null(cache_dir)) {
        cache_dir <- BiocFileCache::getBFCOption("CACHE")
    }
    bfc <- BiocFileCache::BiocFileCache(cache_dir, ask=FALSE)

    url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE132509&format=file"
    gse132509_path <- BiocFileCache::bfcrpath(bfc, "gse132509", fpath=url)
    gse132509_temp <- tempfile()
    untar(gse132509_path, exdir=gse132509_temp)

    print(list.files(gse132509_temp))

    # Load control samples
    sample_files <- list(
        PBMMC_1 = "GSM3872442_PBMMC_1", 
        PBMMC_2 = "GSM3872443_PBMMC_2",
        PBMMC_3 = "GSM3872444_PBMMC_3"
    )
    sobjects <- list()

    for (i in seq_along(sample_files)) {
        sample_name <- names(sample_files)[i]

        mtx_file = paste(sample_files[[i]], "matrix.mtx.gz", sep=".")
        cells_file = paste(sample_files[[i]], "barcodes.tsv.gz", sep=".")
        features_file = paste(sample_files[[i]], "genes.tsv.gz", sep=".")
        counts <- Seurat::ReadMtx(
            mtx = file.path(gse132509_temp, mtx_file),
            cells = file.path(gse132509_temp, cells_file), 
            features = file.path(gse132509_temp, features_file) 
        )

        meta <- data.frame(row.names = colnames(counts))
        meta$sample <- sample_name
        meta$type <- "control"
        sobjects[[sample_name]] <- SeuratObject::CreateSeuratObject(
            counts = counts, min.cells = 3, min.features = 200, meta.data = meta
        )
    }

    # Load case samples
    sample_files <- list(
        PRE_T_1 = "GSM3872440_PRE-T_1", 
        PRE_T_2 = "GSM3872441_PRE-T_2"
    )

    for (i in seq_along(sample_files)) {
        sample_name <- names(sample_files)[i]

        mtx_file = paste(sample_files[[i]], "matrix.mtx.gz", sep=".")
        cells_file = paste(sample_files[[i]], "barcodes.tsv.gz", sep=".")
        features_file = paste(sample_files[[i]], "genes.tsv.gz", sep=".")
        counts <- Seurat::ReadMtx(
            mtx = file.path(gse132509_temp, mtx_file),
            cells = file.path(gse132509_temp, cells_file), 
            features = file.path(gse132509_temp, features_file) 
        )

        meta <- data.frame(row.names = colnames(counts))
        meta$sample <- sample_name
        meta$type <- "case"
        sobjects[[sample_name]] <- SeuratObject::CreateSeuratObject(
            counts = counts, min.cells = 3, min.features = 200, meta.data = meta
        )
    }

    unlink(gse132509_temp, recursive=TRUE)
}