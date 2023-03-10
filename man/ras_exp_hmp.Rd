% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GRN_visualization.R
\name{ras_exp_hmp}
\alias{ras_exp_hmp}
\title{Heatmaps showing mean regulon activity and TF expression of each cluster.}
\usage{
ras_exp_hmp(SeuratObj, rasMat, group.by = "seurat_clusters", assay = "SCT")
}
\arguments{
\item{SeuratObj}{Seurat object}

\item{rasMat}{matrix of regulon activity score(RAS) in each cell}

\item{group.by}{Name of the metadata column to group cells by (for example, seurat_clusters)}

\item{assay}{assay to calculate mean expression}
}
\description{
Heatmaps showing mean regulon activity and TF expression of each cluster.
}
