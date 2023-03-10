% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/integration.R
\name{crossSpecies_integrate}
\alias{crossSpecies_integrate}
\title{Perform integration on datasets of different species with Seurat.}
\usage{
crossSpecies_integrate(
  matrices,
  species,
  resolution = 0.5,
  min.cells = 3,
  min.features = 200
)
}
\arguments{
\item{matrices}{List of single-cell expression matrices of different species.}

\item{species}{species of matrices. Character vector, such as \code{c('Ath', 'Osa')} or \code{c('Ath', 'Osa', 'Zma')},
whose order must match order of matrices. Each word corresponds to each matrix.}

\item{resolution}{resolution parameter for \code{FindClusters}, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities.}

\item{min.cells}{Filter the expression matrices before integration. Filter genes based on the minimum number of cells.}

\item{min.features}{Filter the expression matrices before integration. Filter cells based on minimum number of features.}
}
\value{
a Seurat object
}
\description{
currently integration of three species \code{Arabidopsis thaliana} and \code{Oryza sativa} and \code{Zea mays} are supported.
}
