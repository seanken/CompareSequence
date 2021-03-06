% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Azimuth.R
\name{RunAzimuth_human_cortex}
\alias{RunAzimuth_human_cortex}
\title{Run Azimuth on human cortex}
\usage{
RunAzimuth_human_cortex(seur)
}
\arguments{
\item{seur}{the Seurat object to label cell types for)}
}
\value{
Seurat object with Azimuth results and new assay for SCTransform
}
\description{
Applies azimuth to a Seurat object with the built in human cortex reference
}
