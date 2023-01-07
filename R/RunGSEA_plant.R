
#' Perform Gene Set Enrichment Analysis (GSEA) on Seurat object
#'
#' Perform Gene Set Enrichment Analysis (GSEA) on Seurat object. GSEA is implemented using \code{clusterProfiler} package.
#' This function is only for three plant species: \code{Arabidopsis thaliana} or \code{Oryza sativa} or \code{Zea mays}.
#'
#' @param SeuratObj Seurat object
#' @param by GO KEGG. Will be ignored if parameter "TERM2GENE" is not NULL.
#' @param TERM2GENE Customized terms. Annotation of TERM TO GENE mapping, a data.frame of 2 column with term and gene. NULL by default. If not NULL, parameter "by" will be ignored.
#' @param GeneIDtype Gene ID type in Seurat object. \code{TAIR} for \code{Arabidopsis thaliana}, for example.
#' @param minpct minimum expression percent in each cluster, parameter used to filter Differential expression genes to run GSEA. 0 as default.
#'
#' @importFrom clusterProfiler GSEA gseGO gseKEGG bitr
#' @importFrom magrittr `%>%`
#' @importFrom foreach foreach `%do%`
#'
#' @return Seurat object
#' @export
#'
#' @examples
#' \dontrun{
#' SeuratObj <- RunGSEA_plant(SeuratObj, by = 'GO', GeneIDtype = 'TAIR')
#' }
#'
RunGSEA_plant <- function(SeuratObj,
                          by = 'GO',
                          TERM2GENE = NULL,
                          GeneIDtype = 'TAIR',
                          minpct = 0) {

  by <- match.arg(by, choices = c("GO", "KEGG"))

  if (!"Allmarkers" %in% names(slot(object = SeuratObj, name = 'misc'))) {
    stop("Please run 'RunSeurat()' or 'Seurat::FindAllMarkers()', differentially expressed genes at
         SeuratObj@misc$Allmarkers are required to perform GSEA.")
  }
  Allmarkers <- slot(object = SeuratObj, name = 'misc')[["Allmarkers"]]
  Allmarkers$cluster <- as.character(Allmarkers$cluster)

  if (!is.null(TERM2GENE)) {
    out <- foreach(cls=unique(Allmarkers$cluster)) %do% {
      submarkers <- Allmarkers %>% dplyr::filter(cluster==cls & pct.1 >= minpct)
      geneList <- sort(setNames(submarkers$avg_log2FC, submarkers$gene), decreasing = T)
      res <- GSEA(geneList     = geneList,
                  TERM2GENE    = TERM2GENE,
                  pvalueCutoff = 1,
                  verbose      = FALSE)
      df <- cbind(res@result, cluster=cls)
      df
    }
  } else {

    if (!"species" %in% names(slot(object = SeuratObj, name = 'misc'))) {
      stop("Please run 'readData()' and specify the 'species' parameter.")
    }
    species <- slot(object = SeuratObj, name = 'misc')[["species"]]

    if (is.null(GeneIDtype)) {
      if (!"featureData" %in% names(slot(object = SeuratObj, name = 'misc'))) {
        stop("Please specify the 'GeneIDtype' parameter or run 'DetectGeneIDtype()' to determine the gene ID type of your data.")
      }
      GeneIDtype <- slot(object = SeuratObj, name = 'misc')[["featureData"]][['GeneIDtype']]
    }

    if (by == 'GO') {
      # load the OrgDb
      species <- switch(species,
                        'Arabidopsis thaliana' = 'Ath',
                        'Oryza sativa' = 'Osa',
                        'Zea mays' = 'Zma')
      orgdb <- get_orgDb(species)

      out <- foreach(cls=unique(Allmarkers$cluster)) %do% {
        submarkers <- Allmarkers %>% dplyr::filter(cluster==cls & pct.1 >= minpct)
        geneList <- sort(setNames(submarkers$avg_log2FC, submarkers$gene), decreasing = T)
        res <- gseGO(geneList     = geneList,
                     OrgDb        = orgdb,
                     ont          = "ALL",
                     keyType      = GeneIDtype,
                     nPerm        = 1000,
                     minGSSize    = 10,
                     maxGSSize    = 500,
                     pvalueCutoff = 1,
                     verbose      = FALSE)
        df <- cbind(res@result, cluster=cls)
        df
      }
    }

    if (by == 'KEGG') {
      kegg_code <- CellFunTopic:::get_kegg_code(species)
      # gene ID transition
      if (GeneIDtype != "ENTREZID") {
        species <- switch(species,
                          'Arabidopsis thaliana' = 'Ath',
                          'Oryza sativa' = 'Osa',
                          'Zea mays' = 'Zma')
        orgdb <- get_orgDb(species)
        bitrDF <- bitr(Allmarkers$gene, fromType = GeneIDtype, toType = "ENTREZID", OrgDb = orgdb, drop = TRUE)
        Allmarkers <- merge(Allmarkers, bitrDF, by.x = "gene", by.y = GeneIDtype)
      }

      out <- foreach(cls=unique(Allmarkers$cluster)) %do% {
        submarkers <- Allmarkers %>% dplyr::filter(cluster==cls & pct.1 >= minpct)
        geneList <- sort(setNames(submarkers$avg_log2FC, submarkers$ENTREZID), decreasing = T)
        res <- gseKEGG(geneList          = geneList,
                       organism          = kegg_code,
                       keyType           = 'ncbi-geneid',
                       nPerm             = 1000,
                       minGSSize         = 10,
                       maxGSSize         = 500,
                       pvalueCutoff      = 1,
                       pAdjustMethod     = "BH",
                       verbose           = FALSE,
                       use_internal_data = FALSE)
        df <- cbind(res@result, cluster=cls)
        df
      }
    }
  }

  result <- do.call("rbind", out)
  if (!is.null(TERM2GENE)) {
    slot(object = SeuratObj, name = 'misc')[["GSEAresult_customized"]] <- result
  } else {
    slot(object = SeuratObj, name = 'misc')[[paste0("GSEAresult_", by)]] <- result
  }

  return(SeuratObj)
}










