

#' Toy example single-cell expression matrix of Arabidopsis thaliana.
#'
#' 7000×2000(genes×cells) single-cell expression matrix of Arabidopsis thaliana, taken as a toy example.
#'
#'
#' @format A dgCMatrix of 7000×2000(genes×cells).
#'
"example_Ath"

#' Toy example single-cell expression matrix of Oryza sativa.
#'
#' 7000×2000(genes×cells) single-cell expression matrix of Oryza sativa, taken as a toy example.
#'
#'
#' @format A dgCMatrix of 7000×2000(genes×cells).
#'
"example_Osa"

#' Toy example single-cell expression matrix of Zea mays.
#'
#' 7000×1000(genes×cells) single-cell expression matrix of Zea mays, taken as a toy example.
#'
#'
#' @format A dgCMatrix of 7000×1000(genes×cells).
#'
"example_Zma"


#' One-by-one orthologous genes of rice in Arabidopsis thaliana.
#'
#' One-by-one orthologous genes of rice in Arabidopsis thaliana, used for cross-species scRNA-seq analysis.
#'
#' @keywords internal
#' @format A named vector of 12013 character.
#'
"orthologs_At_Os"


#' One-by-one orthologous genes of Zea mays in Arabidopsis thaliana.
#'
#' One-by-one orthologous genes of Zea mays in Arabidopsis thaliana, used for cross-species scRNA-seq analysis.
#'
#' @keywords internal
#' @format A named vector of 11405 character.
#'
"orthologs_At_Zm"


#' Correspondence between motif and transcription factor in Arabidopsis thaliana.
#'
#' Correspondence between motif and transcription factor in Arabidopsis thaliana.
#'
#' @keywords internal
#' @format A named vector of 1366 character.
#'
"motif2TF_At"


#' Correspondence between motif and transcription factor in Oryza sativa.
#'
#' Correspondence between motif and transcription factor in Oryza sativa.
#'
#' @keywords internal
#' @format A named vector of 476 character.
#'
"motif2TF_Os"


#' Correspondence between motif and transcription factor in Zea mays.
#'
#' Correspondence between motif and transcription factor in Zea mays.
#'
#' @keywords internal
#' @format A named vector of 505 character.
#'
"motif2TF_Zm"


#' GO annotation of Oryza sativa (IRGSP-1.0).
#'
#' A dataset containing rice genes and their GO annotation.
#'
#' @keywords internal
#' @format A data frame with 175728 rows and 3 variables:
#' \describe{
#'   \item{Gene.stable.ID}{gene id of rice}
#'   \item{GO.term.accession}{accession number of GO terms}
#'   \item{GO.term.name}{names of GO terms}
#' }
#' @source \url{https://plants.ensembl.org/index.html}
"GO_rice_IRGSP_1_0"


#' GO annotation of Oryza sativa (v7.0).
#'
#' A dataset containing rice genes and their GO annotation.
#'
#' @keywords internal
#' @format A data frame with 70633  rows and 3 variables:
#' \describe{
#'   \item{Gene.Name}{gene id of rice}
#'   \item{GO.ID}{accession number of GO terms}
#'   \item{GO.Description}{names of GO terms}
#' }
#' @source \url{https://phytozome-next.jgi.doe.gov/}
"GO_rice_V7"


#' GO annotation of Zea mays (Zm-B73-REFERENCE-NAM-5.0).
#'
#' A dataset containing Zea mays genes and their GO annotation.
#'
#' @keywords internal
#' @format A data frame with 291353 rows and 3 variables:
#' \describe{
#'   \item{Gene.stable.ID}{gene id of Zea mays}
#'   \item{GO.term.accession}{accession number of GO terms}
#'   \item{GO.term.name}{names of GO terms}
#' }
#' @source \url{https://plants.ensembl.org/index.html}
"GO_zeamays_NAM_5"


#' GO annotation of Zea mays (RefGen_V4).
#'
#' A dataset containing Zea mays genes and their GO annotation.
#'
#' @keywords internal
#' @format A data frame with 233821 rows and 3 variables:
#' \describe{
#'   \item{Gene.Name}{gene id of Zea mays}
#'   \item{GO.ID}{accession number of GO terms}
#'   \item{GO.Description}{names of GO terms}
#' }
#' @source \url{https://phytozome-next.jgi.doe.gov/}
"GO_zeamays_V4"


#' KEGG annotation of Oryza sativa (v7.0).
#'
#' A dataset containing rice genes and their KEGG annotation.
#'
#' @keywords internal
#' @format A data frame with 10186  rows and 3 variables:
#' \describe{
#'   \item{Gene.Name}{gene id of rice}
#'   \item{KEGG.ID}{accession number of KEGG terms}
#'   \item{KEGG.Description}{names of KEGG terms}
#' }
#' @source \url{https://phytozome-next.jgi.doe.gov/}
"KEGG_rice_V7"

#' KEGG annotation of Zea mays (Zm-B73-REFERENCE-NAM-5.0).
#'
#' A dataset containing Zea mays genes and their KEGG annotation.
#'
#' @keywords internal
#' @format A data frame with 16073  rows and 3 variables:
#' \describe{
#'   \item{Gene.Name}{gene id of Zea mays}
#'   \item{KEGG.ID}{accession number of KEGG terms}
#'   \item{KEGG.Description}{names of KEGG terms}
#' }
#' @source \url{https://phytozome-next.jgi.doe.gov/}
"KEGG_zeamays_NAM_5"

#' KEGG annotation of Zea mays (RefGen_V4).
#'
#' A dataset containing Zea mays genes and their KEGG annotation.
#'
#' @keywords internal
#' @format A data frame with 22136 rows and 3 variables:
#' \describe{
#'   \item{Gene.Name}{gene id of Zea mays}
#'   \item{KEGG.ID}{accession number of KEGG terms}
#'   \item{KEGG.Description}{names of KEGG terms}
#' }
#' @source \url{https://phytozome-next.jgi.doe.gov/}
"KEGG_zeamays_V4"

















