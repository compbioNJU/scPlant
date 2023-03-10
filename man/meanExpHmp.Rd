% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/decomposeBulk.R
\name{meanExpHmp}
\alias{meanExpHmp}
\title{Heatmap showing mean expression of top 20 markers of each cell type in each sample.}
\usage{
meanExpHmp(SeuratObj, exprs, meta)
}
\arguments{
\item{SeuratObj}{A Seurat object.}

\item{exprs}{A data frame or matrix of raw read counts of bulk RNA-seq samples. Column names correspond to sample names, row names to genes.}

\item{meta}{A data frame of 2 columns with Sample and Class, correlate sample names with sample class.}
}
\description{
Heatmap showing mean expression of top 20 markers of each cell type in each sample.
}
