% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Auto_annotation.R
\name{AutoAnnotate_Garnett}
\alias{AutoAnnotate_Garnett}
\title{Automatic annotation of Seurat object using Garnett package.}
\usage{
AutoAnnotate_Garnett(
  SeuratObj,
  marker_file_path,
  gene_id_type,
  species = c("Ath", "Osa", "Zma")
)
}
\arguments{
\item{SeuratObj}{Seurat object to annotate}

\item{marker_file_path}{see https://cole-trapnell-lab.github.io/garnett/docs/ to construct a marker file}

\item{gene_id_type}{gene ID type of Seurat object and marker file}

\item{species}{\code{c('Ath', 'Osa', 'Zma')}. Currently three plant species are supported.}
}
\value{
Seurat object. The predicted cell type label is stored in \code{SeuratObj$predicted_label}.
}
\description{
Automatic annotation of Seurat object using Garnett package.
}
