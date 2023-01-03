
#' Draw Constellation plot with scrattch.hicat package.
#'
#' @param SeuratObj Seurat object
#' @param cluster_label one column of \code{SeuratObj@meta.data}, used for cluster annotation to display in the plot.
#' @param out.dir location to write plotting files to
#' @param k K for KNN algorithm.
#' @param knn.outlier.th  Threshold to determine if a nearest neighbor is too far
#' @param outlier.frac.th Threshold to determine if a cell is a outlier if it is too far away from most of its neighbors.
#'
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' SeuratObj <- pbmc3k.SeuratData::pbmc3k.final
#' Seurat::Idents(SeuratObj) <- SeuratObj$seurat_clusters
#' constellationPlot(SeuratObj, cluster_label = "seurat_annotations", out.dir = "Constellation_plots")
#' }
#'
constellationPlot <- function(SeuratObj, cluster_label = "cell_type", out.dir = "Constellation_plots",
                              k=15, knn.outlier.th=2, outlier.frac.th=0.5) {
  # get output of KNN.graph
  rd.dat <- Seurat::Embeddings(SeuratObj, reduction = 'pca')
  cl <- setNames(as.character(Seurat::Idents(SeuratObj)), names(Seurat::Idents(SeuratObj)))
  cl.df <- unique(data.frame(cluster=cl,
                             cluster_label=as.character(SeuratObj@meta.data[names(Seurat::Idents(SeuratObj)), cluster_label]),
                             stringsAsFactors = F))
  rownames(cl.df) <- cl.df$cluster
  knnResult <- scrattch.hicat::get_knn_graph(rd.dat, cl, cl.df, k=k, knn.outlier.th=knn.outlier.th, outlier.frac.th=outlier.frac.th)
  knn.cl.df <- knnResult[['knn.cl.df']]
  knn.cl.df$cl.from <- as.character(knn.cl.df$cl.from)
  knn.cl.df$cl.to <- as.character(knn.cl.df$cl.to)
  # Draw Constellation plot
  cell_eb <- as.data.frame(Seurat::Embeddings(SeuratObj, reduction = 'umap')[,1:2])
  cell_eb$cl <- as.character(Seurat::Idents(SeuratObj)[rownames(cell_eb)])
  colnames(cell_eb) <- c('x','y','cl')
  cl.center.df <- cell_eb %>% group_by(cl) %>% dplyr::summarise(x=median(x), y=median(y), .groups = "keep") %>% as.data.frame
  cl.center.df$cluster_id <- cl.center.df$cl
  cl.center.df$cluster_color <- CellFunTopic::scPalette2(length(unique(cl.center.df$cl)))
  nn <- tibble::deframe(cl.df)
  cl.center.df$cluster_label <- nn[as.character(cl.center.df$cl)]
  cl.center.df$cluster_size <- as.numeric(table(cl)[as.character(cl.center.df$cl)])

  plotting.constellation <- scrattch.hicat::plot_constellation(knn.cl.df = knn.cl.df, cl.center.df = cl.center.df,
                                                               out.dir = out.dir, node.dodge=TRUE)
}








