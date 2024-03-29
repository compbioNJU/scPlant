% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MINIEX_visualization.R
\name{targets_MINIEX}
\alias{targets_MINIEX}
\title{Network diagram showing targets of each regulon}
\usage{
targets_MINIEX(
  MINIEXouputPath,
  cluster = NULL,
  regulons = NULL,
  tf.size = 7,
  target.size = 4,
  tf.label.cex = 1,
  target.label.cex = 0.5
)
}
\arguments{
\item{MINIEXouputPath}{output path of MINI-EX.}

\item{cluster}{which cluster to show.}

\item{regulons}{regulons to show.}

\item{tf.size}{size of TF nodes}

\item{target.size}{size of target nodes}

\item{tf.label.cex}{label size of TF nodes}

\item{target.label.cex}{label size of target nodes}
}
\description{
Network diagram showing targets of each regulon inferred by MINI-EX.
}
\examples{
\dontrun{
targets_MINIEX(MINIEXouputPath="regulons_output", cluster = 1,
               regulons = c("AT4G25490", "AT3G58710"))
}


}
