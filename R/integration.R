

#' Perform integration on datasets of different species with Seurat.
#'
#' currently integration of three species \code{Arabidopsis thaliana} and \code{Oryza sativa} and \code{Zea mays} are supported.
#'
#' @param matrices List of single-cell expression matrices of different species.
#' @param species species of matrices. Character vector, such as \code{c('Ath', 'Osa')} or \code{c('Ath', 'Osa', 'Zma')},
#' whose order must match order of matrices. Each word corresponds to each matrix.
#' @param resolution resolution parameter for \code{FindClusters}, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities.
#' @param min.cells Filter the expression matrices before integration. Filter genes based on the minimum number of cells.
#' @param min.features Filter the expression matrices before integration. Filter cells based on minimum number of features.
#'
#' @importFrom Seurat CreateSeuratObject SCTransform SelectIntegrationFeatures PrepSCTIntegration FindIntegrationAnchors
#' @importFrom Seurat IntegrateData RunPCA RunUMAP FindNeighbors FindClusters
#'
#' @return a Seurat object
#' @export
#'
crossSpecies_integrate <- function(matrices, species, resolution = 0.5, min.cells = 3, min.features = 200) {

  # filter the expression matrices
  matrices <- lapply(matrices, filterMat, min.cells = min.cells, min.features = min.features)
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
  genes <- Reduce(intersect, lapply(obj_list, rownames))
  obj_list <- lapply(obj_list, function(x){
    subset(x, features = genes)
  })
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
  # Filter cells based on minimum number of features
  if (min.features > 0) {
    nfeatures <- Matrix::colSums(x = counts > 0)
    counts <- counts[, which(x = nfeatures >= min.features)]
  }
  # Filter genes based on the minimum number of cells
  if (min.cells > 0) {
    num.cells <- Matrix::rowSums(x = counts > 0)
    counts <- counts[which(x = num.cells >= min.cells), ]
  }
  return(counts)
}


#' Bar plot showing the percentage of cells from different species.
#'
#' Bar plot showing the percentage of cells from different species after cross-species integration.
#'
#' @param SeuratObj Seurat object.
#' @param group_by column of \code{SeuratObj@meta.data}.
#'
#' @importFrom ggplot2 ggplot aes geom_bar theme_classic xlab theme element_text
#'
#' @return a ggplot ibject.
#' @export
#'
species_percentage <- function(SeuratObj, group_by = 'seurat_clusters') {
  mm <- table(as.character(SeuratObj@meta.data[[group_by]]), SeuratObj$species)
  pp <- as.data.frame(mm) %>% magrittr::set_colnames(c('group', 'species', 'cellnumber')) %>%
    ggplot(aes(x=cellnumber, y=group, fill=species)) + geom_bar(stat = 'identity', position = 'fill') +
    theme_classic() + xlab('Percentage(%)') + theme(axis.text = element_text(color='black'))
  return(pp)
}

















