

#' Perform integration on datasets of different species with Seurat.
#'
#' currently integration of three species \code{Arabidopsis thaliana} and \code{Oryza sativa} and \code{Zea mays} are supported.
#'
#' @param matrices List of single-cell expression matrices of different species.
#' @param species species of matrices. Charactor vector, such as \code{c('Ath', 'Osa')} or \code{c('Ath', 'Osa', 'Zma')},
#' whose order must match order of matrices..
#' @param resolution resolution parameter for \code{FindClusters}, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities.
#'
#' @importFrom Seurat CreateSeuratObject SCTransform SelectIntegrationFeatures PrepSCTIntegration FindIntegrationAnchors
#' @importFrom Seurat IntegrateData RunPCA RunUMAP FindNeighbors FindClusters
#'
#' @return
#' @export
#'
crossSpecies_integrate <- function(matrices, species, resolution = 0.5) {

  # filter the expression matrices
  matrices <- lapply(matrices, filterMat)
  # transform genes to orthologous gene in Arabidopsis thaliana and create Seurat object
  athlist <- osalist <- zmalist <- NULL
  if ('Ath' %in% species) {
    athlist <- lapply(matrices[species == 'Ath'], function(mm){
      colnames(mm) <- paste0("ath@", colnames(mm))
      object <- CreateSeuratObject(counts = mm, project = "Ath", min.cells = 0, min.features = 0)
      object@meta.data$species <- 'Arabidopsis thaliana'
      object
    })
  }
  if ('Osa' %in% species) {
    osalist <- lapply(matrices[species == 'Osa'], function(mm){
      colnames(mm) <- paste0("osa@", colnames(mm))
      mm <- mm[rownames(mm) %in% unique(names(orthologs_At_Os)), ]
      rownames(mm) <- unname(orthologs_At_Os[rownames(mm)])
      object <- CreateSeuratObject(counts = mm, project = "Osa", min.cells = 0, min.features = 0)
      object@meta.data$species <- 'Oryza sativa'
      object
    })
  }
  if ('Zma' %in% species) {
    zmalist <- lapply(matrices[species == 'Zma'], function(mm){
      colnames(mm) <- paste0("zma@", colnames(mm))
      mm <- mm[rownames(mm) %in% unique(names(orthologs_At_Zm)), ]
      rownames(mm) <- unname(orthologs_At_Zm[rownames(mm)])
      object <- CreateSeuratObject(counts = mm, project = "Zma", min.cells = 0, min.features = 0)
      object@meta.data$species <- 'Zea mays'
      object
    })
  }
  obj_list <- c(athlist, osalist, zmalist)
  # Performing integration on datasets normalized with SCTransform
  obj_list <- lapply(X = obj_list, FUN = SCTransform)
  features <- SelectIntegrationFeatures(object.list = obj_list, nfeatures = 1000)
  obj_list <- PrepSCTIntegration(object.list = obj_list, anchor.features = features)
  anchors <- FindIntegrationAnchors(object.list = obj_list, normalization.method = "SCT",
                                           anchor.features = features)
  combined.sct <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
  combined.sct <- RunPCA(combined.sct, verbose = FALSE)
  combined.sct <- RunUMAP(combined.sct, reduction = "pca", dims = 1:30)
  combined.sct <- FindNeighbors(combined.sct, reduction = "pca", dims = 1:30)
  combined.sct <- FindClusters(combined.sct, resolution = resolution)
  return(combined.sct)
}




filterMat <- function(counts, min.cells = 3, min.features = 200) {
  # Filter based on min.features
  if (min.features > 0) {
    nfeatures <- Matrix::colSums(x = counts > 0)
    counts <- counts[, which(x = nfeatures >= min.features)]
  }
  # filter genes on the number of cells expressing
  if (min.cells > 0) {
    num.cells <- Matrix::rowSums(x = counts > 0)
    counts <- counts[which(x = num.cells >= min.cells), ]
  }
  return(counts)
}

