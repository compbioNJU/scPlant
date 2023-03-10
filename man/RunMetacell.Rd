% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RunMetacell.R
\name{RunMetacell}
\alias{RunMetacell}
\title{Run metacell on Seurat object}
\usage{
RunMetacell(SeuratObj, species)
}
\arguments{
\item{SeuratObj}{Seurat object}

\item{species}{species}
}
\description{
The MetaCell R package facilitates analysis of single cell RNA-seq UMI matrices by computing partitions
of a cell similarity graph into small (~20-200 typically) homogeneous groups of cells which are defined as
metacells (MCs). The derived MCs are then used for building different representations of the data,
allowing matrix or 2D graph visualization forming a basis for analysis of cell types, subtypes,
transcriptional gradients, cell-cycle variation, gene modules and their regulatory models and more.
More details on the usage of the MetaCell pipeline is available in the package vignettes \url{https://tanaylab.github.io/metacell/index.html}.
Note that metacell package is tested on linux and macbooks, and is currently not compatible on Windows.
}
