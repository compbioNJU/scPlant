% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GRN_visualization.R
\name{topRegulons}
\alias{topRegulons}
\title{Network diagram showing top regulons of each cluster}
\usage{
topRegulons(
  rssMat,
  topn = 5,
  from.size = 7,
  to.size = 4,
  from.label.cex = 1,
  to.label.cex = 0.7
)
}
\arguments{
\item{rssMat}{matrix of regulon specificity score(RSS) in each cluster}

\item{topn}{number of top regulons to draw}

\item{from.size}{size of 'cluster' nodes}

\item{to.size}{size of 'regulon' nodes}

\item{from.label.cex}{label size of 'cluster' nodes}

\item{to.label.cex}{label size of 'regulon' nodes}
}
\description{
Network diagram showing top regulons of each cluster according to regulon specificity score(RSS).
}
