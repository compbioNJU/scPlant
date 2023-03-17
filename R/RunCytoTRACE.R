
#' Run CytoTRACE on Seurat object
#'
#' @param SeuratObj A Seurat object
#' @param keep_orig_reduction  keep original reduction dimension of Seurat object.
#' @param reduction_method reduction method, umap as default. Only works when \code{keep_orig_reduction = TRUE}.
#' @param ncores the number of cores to utilize.
#' @param save_plot  save plots to pdf or not.
#' @param return_object return CytoTRACE result or not.
#'
#' @return CytoTRACE result. See \code{?CytoTRACE::CytoTRACE}.
#' @export
#'
RunCytoTRACE <- function(SeuratObj, 
                         keep_orig_reduction = TRUE,
                         reduction_method = "umap",
                         ncores = 1,
                         save_plot = TRUE,
                         return_object = TRUE) {
  
  mm <- as.matrix(SeuratObj@assays$RNA@counts)
  results <- CytoTRACE::CytoTRACE(mat = mm, ncores = ncores, subsamplesize = 1000)
  
  if (keep_orig_reduction) {
    if (save_plot) {
      if (!dir.exists(paths = "./CytoTRACEoutput")) {
        dir.create("./CytoTRACEoutput", recursive = TRUE)
      }
      saveRDS(results, './CytoTRACEoutput/cytotrace_result.rds')
      CytoTRACE::plotCytoTRACE(results, outputDir = "./CytoTRACEoutput/", emb = Seurat::Embeddings(SeuratObj, reduction=reduction_method))
      CytoTRACE::plotCytoGenes(results, numOfGenes = 10, outputDir = "./CytoTRACEoutput/")
    }
    CytoTRACE::plotCytoTRACE(results, emb = Seurat::Embeddings(SeuratObj, reduction=reduction_method))
  } else {
    if (save_plot) {
      if (!dir.exists(paths = "./CytoTRACEoutput")) {
        dir.create("./CytoTRACEoutput", recursive = TRUE)
      }
      saveRDS(results, './CytoTRACEoutput/cytotrace_result.rds')
      CytoTRACE::plotCytoTRACE(results, outputDir = "./CytoTRACEoutput/")
      CytoTRACE::plotCytoGenes(results, numOfGenes = 10, outputDir = "./CytoTRACEoutput/")
    }
    CytoTRACE::plotCytoTRACE(results)
  }
  
  if (return_object) {
    return(results)
  }
}













