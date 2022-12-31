


#' Run slingshot on Seurat object
#'
#' @param SeuratObj A Seurat object
#' @param clusterlabel cluster labels
#' @param reduction reduction
#' @param dynamic_genes whether to identify dynamic genes across pseudotime
#' @param save_plot save plots to pdf or not.
#' @param return_object return object or not.
#'
#' @importFrom slingshot slingshot SlingshotDataSet
#' @importFrom SingleCellExperiment reducedDims
#' @importFrom stats heatmap
#'
#' @return if parameter \code{return_object = TRUE}, SingleCellExperiment object and dynamic genes will be returned.
#' @export
#'
RunSlingshot <- function(SeuratObj,
                         clusterlabel = 'seurat_clusters',
                         reduction = 'UMAP',
                         dynamic_genes = FALSE,
                         save_plot = TRUE,
                         return_object = TRUE) {

  sce <- Seurat::as.SingleCellExperiment(SeuratObj)
  sce <- slingshot(sce, clusterLabels = clusterlabel, reducedDim = reduction)

  if (save_plot) {
    if (!dir.exists(paths = "./output/plots/RunSlingshot")) {
      dir.create("./output/plots/RunSlingshot", recursive = TRUE)
    }
    pdf(file = "./output/plots/RunSlingshot/RunSlingshot_1.pdf", width = 10, height = 10)
    colors <- grDevices::colorRampPalette(pals::brewer.spectral(11)[-6])(100)
    plotcol <- colors[cut(sce$slingPseudotime_1, breaks=100)]
    plot(reducedDims(sce)[[reduction]], col = plotcol, pch=16, asp = 1)
    lines(SlingshotDataSet(sce), lwd=2, col='black')
    dev.off()
    pdf(file = "./output/plots/RunSlingshot/RunSlingshot_2.pdf", width = 10, height = 10)
    clustercolors <- setNames(CellFunTopic::scPalette2(length(unique(sce[[clusterlabel]]))), unique(sce[[clusterlabel]]))
    plot(reducedDims(sce)[[reduction]], col = clustercolors[as.character(sce[[clusterlabel]])], pch=16, asp = 1)
    lines(SlingshotDataSet(sce), lwd=2, type = 'lineages', col = 'black')
    dev.off()
  }

  if (dynamic_genes) {
    # fit negative binomial GAM
    sce <- tradeSeq::fitGAM(sce)
    # test for dynamic expression
    ATres <- tradeSeq::associationTest(sce)
    if (save_plot) {
      # pick out the top 250 genes based on p-values and visualize their expression over developmental time with a heatmap.
      pdf(file = "./output/plots/RunSlingshot/dynamic_genes.pdf", width = 10, height = 10)
      topgenes <- rownames(ATres[order(ATres$pvalue), ])[1:250]
      pst.ord <- order(sce$slingPseudotime_1, na.last = NA)
      heatdata <- assays(sce)$counts[topgenes, pst.ord]
      heatclus <- sce[[clusterlabel]][pst.ord]
      heatmap(log1p(heatdata), Colv = NA, ColSideColors = clustercolors[as.character(heatclus)])
      dev.off()
    }
  }
  if (return_object) {
    if (dynamic_genes) {
      return(list(SCE = sce, dynamic_genes = ATres))
    } else {
      return(list(SCE = sce))
    }
  }
}






