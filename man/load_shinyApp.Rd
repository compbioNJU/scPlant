% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/load_shinyApp.R
\name{load_shinyApp}
\alias{load_shinyApp}
\title{visualize the analysis result in a built-in shiny app}
\usage{
load_shinyApp(SeuratObj, rasMat, rssMat, tf_target)
}
\arguments{
\item{SeuratObj}{Seurat object}

\item{rasMat}{matrix of regulon activity score(RAS) in each cell}

\item{rssMat}{matrix of regulon specificity score(RSS) in each cluster}

\item{tf_target}{TF and targets of each regulon. This information comes from post-processing of pyscenic.}
}
\description{
visualize the analysis result in a built-in shiny app
}
\examples{
\dontrun{
if(interactive()) {
    load_shinyApp(SeuratObj, rasMat, rssMat, tf_target)
}
}

}
