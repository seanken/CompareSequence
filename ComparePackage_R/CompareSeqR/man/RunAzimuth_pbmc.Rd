% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Azimuth.R
\name{RunAzimuth_pbmc}
\alias{RunAzimuth_pbmc}
\title{Run Azimuth on PBMC}
\usage{
RunAzimuth_pbmc(seur, col = "celltype.l2")
}
\arguments{
\item{seur}{the Seurat object}

\item{col}{the column in the meta data to use (celltype.l1, celltype.l2, or celltype.l3 all being difference options)}
}
\value{
Seurat object with Azimuth results and new assay for SCTransform
}
\description{
Applies azimuth to a Seurat object with the built in human PBMC reference
}
