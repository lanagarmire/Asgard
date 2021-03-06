% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SCalignment.R
\name{SCalignment}
\alias{SCalignment}
\title{Single-cell Alignment.}
\usage{
SCalignment(
  SC.list = SC.list,
  CellCycle = TRUE,
  anchor.features = 2000,
  by.CellType = FALSE
)
}
\arguments{
\item{SC.list}{A list contains single-cell data.}

\item{CellCycle}{logical; if TRUE, will do cell-cycle regression. The default value is TURE.}

\item{anchor.features}{The number of "anchor" genes used in single-cell pairwise correspondence. The default value is 2000.}

\item{by.CellType}{logical; if TRUE, will annotate cell type. The default value is FALSE.}
}
\value{
A Seurat object of aligned cells.
}
\description{
Single-cell pairwise correspondence, clustering and annotation.
}
\details{
In Asgard, every cell type in the diseased sample is paired to that in the normal (or control) sample, according to “anchor” genes that are consistently expressed between diseased and normal cells. It then clusters paired cells and identifies differentially expressed genes between paired diseased and normal cells for every clustered single cell population. Alternatively, in a cell-type-based procedure, Asgard annotates cell type for every single cell population and identifies differentially expressed genes between paired diseased and normal cells for every type of cell.
}
