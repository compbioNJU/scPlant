#' Network diagram showing targets of each regulon
#'
#' Network diagram showing targets of each regulon inferred by MINI-EX.
#'
#' @param MINIEXouputPath output path of MINI-EX.
#' @param cluster which cluster to show.
#' @param regulons regulons to show.
#' @param tf.size size of TF nodes
#' @param target.size size of target nodes
#' @param tf.label.cex label size of TF nodes
#' @param target.label.cex label size of target nodes
#'
#' @importFrom igraph graph_from_data_frame E V `E<-` `V<-`
#'
#' @export
#'
#' @examples
#' \dontrun{
#' targets_MINIEX(MINIEXouputPath="regulons_output", cluster = 1,
#'                regulons = c("AT4G25490", "AT3G58710"))
#' }
#'
#'
targets_MINIEX <- function(MINIEXouputPath, cluster = NULL, regulons = NULL,
                           tf.size = 7, target.size = 4,
                           tf.label.cex = 1, target.label.cex = 0.5) {

  if (length(cluster)>1) {
    stop("Only show regulons of one cluster at once. Please offer one cluster at each run.")
  }
  output <- read.table(file = list.files(path = MINIEXouputPath, pattern = '_regulons.txt', full.names = T), stringsAsFactors = F)
  cluster <- paste0('Cluster_', cluster)
  tf_target <- output %>% dplyr::filter(V1 %in% regulons & V2 == cluster)
  if (nrow(tf_target) < 1) {
    stop("No result with the parameter regulons/cluster provided. Please check your parameter settings.")
  }
  ll <- lapply(1:nrow(tf_target), function(x){
    targets <- tf_target[x, 3]
    targets <- strsplit(targets, split = ",")[[1]]
    data.frame(source = tf_target[x, 1], target=targets)
  })
  edges <- do.call('rbind', ll)

  net <- graph_from_data_frame(d=edges, directed=T)
  E(net)$width <- 0.2
  E(net)$arrow.size <- .1
  E(net)$arrow.width <- .1
  V(net)$label.cex <- ifelse(V(net)$name %in% unique(edges$source), tf.label.cex, target.label.cex)
  V(net)$label.font <- 4
  V(net)$label.color <- "black"
  V(net)$frame.color <- NA
  V(net)$color <- ifelse(V(net)$name %in% unique(edges$source), pals::brewer.blues(9)[4], pals::brewer.greys(9)[3])
  V(net)$size <- ifelse(V(net)$name %in% unique(edges$source), tf.size, target.size)
  V(net)$shape <- 'circle'
  set.seed(123)
  plot(net, layout=igraph::layout_with_kk, vertex.label.degree=0, vertex.label.dist=0)
}



#' Network diagram showing top regulons of each cluster
#'
#' Network diagram showing top regulons of each cluster according to Borda ranking.
#'
#' @param MINIEXouputPath output path of MINI-EX.
#' @param topn number of top regulons to draw
#' @param from.size size of 'cluster' nodes
#' @param to.size size of 'regulon' nodes
#' @param from.label.cex label size of 'cluster' nodes
#' @param to.label.cex label size of 'regulon' nodes
#'
#'
#' @importFrom igraph graph_from_data_frame E V `E<-` `V<-`
#' @importFrom stats setNames
#'
#' @export
#'
#' @examples
#' \dontrun{
#' topRegulons_MINIEX(MINIEXouputPath="regulons_output", topn = 5)
#' }
#'
topRegulons_MINIEX <- function(MINIEXouputPath, topn = 5,
                               from.size = 7, to.size = 3,
                               from.label.cex = 1, to.label.cex = 0.5) {

  # choose top regulons of each cluster according to Borda ranking
  output <- readxl::read_xlsx(list.files(path = MINIEXouputPath, pattern = '_rankedRegulons.xlsx', full.names = T))
  output <- output[, c("cluster", "TF", "borda_clusterRank")]
  chosen_regulon <- output %>% dplyr::group_by(cluster) %>%
    dplyr::slice_min(order_by = borda_clusterRank, n = topn, with_ties = F) %>% dplyr::ungroup()

  # network diagram
  edgedf <- chosen_regulon[, c("cluster", "TF")] %>% magrittr::set_colnames(c('source', 'target'))
  celltypes <- unique(edgedf$source)
  nodeColor <- setNames(CellFunTopic::scPalette2(length(celltypes)), celltypes)
  edgedf$color = unname(nodeColor[edgedf$source])

  net <- graph_from_data_frame(edgedf, directed = T)
  E(net)$arrow.size <- .2
  E(net)$arrow.width <- .2
  V(net)$label.cex <- ifelse(V(net)$name %in% celltypes, from.label.cex, to.label.cex)
  V(net)$label.font <- ifelse(V(net)$name %in% celltypes, 1, 4)
  V(net)$label.color <- "black"
  V(net)$frame.color <- NA
  V(net)$size <- ifelse(V(net)$name %in% celltypes, from.size, to.size)
  V(net)$shape <- "circle"
  V(net)$color <- ifelse(V(net)$name %in% celltypes, nodeColor[V(net)$name], pals::brewer.greys(9)[3])
  V(net)$label.dist <- ifelse(V(net)$name %in% celltypes, 0, 0.5)
  set.seed(123)
  plot(net, vertex.label.degree=0)
}



#' Dot plot showing top regulons of each cluster according to Borda ranking
#'
#' @param MINIEXouputPath output path of MINI-EX.
#' @param cluster which cluster to show
#' @param topn number of top regulons to draw
#'
#' @importFrom ggplot2 geom_point ylab xlab theme_bw theme element_blank element_line element_text
#'
#' @export
#'
#' @examples
#' \dontrun{
#' BordaRank_MINIEX(MINIEXouputPath="regulons_output", cluster=1, topn = 10)
#' }
#'
BordaRank_MINIEX <- function(MINIEXouputPath, cluster, topn = 10) {

  # choose top regulons of each cluster according to Borda ranking
  output <- readxl::read_xlsx(list.files(path = MINIEXouputPath, pattern = '_rankedRegulons.xlsx', full.names = T))
  output <- output[, c("cluster", "TF", "borda_clusterRank")]
  cluster <- paste0('Cluster_', cluster)

  data <- output[stringr::str_ends(string = output$cluster, pattern = cluster), ] %>%
    dplyr::select(TF, borda_clusterRank) %>%
    dplyr::arrange(borda_clusterRank) %>% tibble::rowid_to_column("index")

  data$pt.col <- ifelse(data$index <= topn, "#007D9B", "#BECEE3")
  data <- head(data, n=200)
  data.label <- head(data, n=topn)

  pp <- ggplot(data, aes(index, borda_clusterRank)) +
    geom_point(size=3, color=data$pt.col) +
    ggrepel::geom_text_repel(inherit.aes = FALSE, data = data.label, aes(index, borda_clusterRank, label=TF), size=4,
                             max.overlaps = 30) +
    ggtitle(as.character(cluster)) + ylab("Borda ranking") + xlab('Regulons') +
    theme_bw(base_size = 12) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(color="black"),
          axis.ticks = element_line(color="black"),
          axis.text = element_text(color = "black"),
          plot.title = element_text(hjust = .5)
    )
  return(pp)
}


#' Heatmaps showing cluster enrichment and TF expression of each regulon in each cluster.
#'
#' @param SeuratObj Seurat object
#' @param MINIEXouputPath output path of MINI-EX.
#' @param group.by Name of the metadata column to group cells by (for example, seurat_clusters)
#' @param assay assay to calculate mean expression
#'
#' @importFrom ComplexHeatmap Heatmap draw
#' @importFrom grid unit
#' @importFrom dplyr `%>%`
#'
#' @export
#'
#' @examples
#' \dontrun{
#' enrich_exp_hmp(SeuratObj, MINIEXouputPath="regulons_output")
#' }
#'
enrich_exp_hmp <- function(SeuratObj, MINIEXouputPath, group.by = "seurat_clusters", assay = 'SCT') {

  # cluster enrichment of each regulon
  output <- readxl::read_xlsx(list.files(path = MINIEXouputPath, pattern = '_rankedRegulons.xlsx', full.names = T))
  output <- output[, c("cluster", "TF", "qval_cluster")]
  output$enrichment <- -log10(output$qval_cluster)
  enrichment <- reshape2::acast(output, TF~cluster, value.var = 'enrichment')
  enrichment[is.na(enrichment)] <- 0
  enrichment <- t(scale(t(enrichment), center = T, scale=T))
  col_fun1 = circlize::colorRamp2(seq(min(enrichment), max(enrichment), length.out = 10), c("white", pals::brewer.orrd(9)))
  ph <- pheatmap::pheatmap(enrichment, cluster_rows = T, cluster_cols = T, silent = T)
  rowOrder <- ph$tree_row$labels[ph$tree_row$order]
  colOrder <- ph$tree_col$labels[ph$tree_col$order]

  # mean TF expression
  Seurat::Idents(SeuratObj) <- SeuratObj@meta.data[[group.by]]
  meanExp <- Seurat::AverageExpression(SeuratObj, assays = assay, slot = 'data')[[assay]]
  meanExp <- meanExp[rownames(enrichment), ] %>% as.matrix()
  meanExp <- t(scale(t(meanExp), center = T, scale=T))
  col_fun2 = circlize::colorRamp2(seq(min(meanExp), max(meanExp), length.out = 10), c("white", pals::brewer.blues(9)))

  # column names(cluster names) may be different.
  if (all(colnames(enrichment) %in% colnames(meanExp))) {
    oo <- colOrder
  } else {
    colnames(meanExp) <- paste0('Cluster_', colnames(meanExp))
    oo <- stringr::str_extract(colOrder, pattern = "Cluster_.+")
  }

  ht1 <- Heatmap(enrichment[rowOrder, colOrder], name = "Enrichment", col = col_fun1,
                 column_names_rot = 45, row_names_gp = gpar(fontsize = 1),
                 cluster_rows = F, cluster_columns = F,
                 column_names_gp = gpar(fontsize = 7), column_names_centered = F,  border = TRUE)
  ht2 <- Heatmap(meanExp[rowOrder, oo], name = "TF Expression", col = col_fun2,
                 column_names_rot = 45, row_names_gp = gpar(fontsize = 1),
                 cluster_rows = F, cluster_columns = F,
                 column_names_gp = gpar(fontsize = 7), column_names_centered = F,  border = TRUE)
  ht_list = ht1 + ht2
  draw(ht_list, ht_gap = unit(1, "cm"), heatmap_legend_side = "left")
}



#' Dimension reduction plot showing cluster enrichment and TF expression
#'
#' @param SeuratObj Seurat object
#' @param MINIEXouputPath output path of MINI-EX.
#' @param gene TF to show
#' @param reduction umap or tsne
#'
#' @importFrom Seurat FeaturePlot NoAxes
#' @importFrom ggplot2 ggtitle guides guide_colorbar
#'
#' @export
#'
#' @examples
#' \dontrun{
#' enrich_exp_scatter(SeuratObj, MINIEXouputPath="regulons_output", gene = 'AT2G40750')
#' }
#'
enrich_exp_scatter <- function(SeuratObj, MINIEXouputPath, gene, reduction = 'umap') {

  p1 <- Seurat::FeaturePlot(SeuratObj, reduction = reduction, order=TRUE,
                            cols=c('lightgrey', pals::brewer.blues(10)),
                            coord.fixed = F,
                            pt.size=0.2, combine=T, features = gene) + Seurat::NoAxes() + #NoLegend() +
    # theme(plot.title = element_blank())
    ggtitle(gene) + guides(color=guide_colorbar(title="Expression"))

  # cluster enrichment of each regulon
  output <- readxl::read_xlsx(list.files(path = MINIEXouputPath, pattern = '_rankedRegulons.xlsx', full.names = T))
  output <- output[, c("cluster", "TF", "qval_cluster")]
  output$enrichment <- -log10(output$qval_cluster)
  output$cluster <- stringr::str_extract(output$cluster, pattern = "Cluster_.+") %>%
    stringr::str_split(pattern = "_", simplify = T) %>% .[, 2]
  nn <- output %>% dplyr::filter(TF == gene) %>% dplyr::select(cluster, enrichment) %>%
    unique %>% tibble::deframe()
  SeuratObj$Regulon <- unname(nn[as.character(SeuratObj$seurat_clusters)])
  p2 <- Seurat::FeaturePlot(SeuratObj, reduction = reduction, order=TRUE, coord.fixed = F,
                            cols=c("lightgrey", pals::brewer.orrd(10)),
                            pt.size=0.2, combine=T,
                            features="Regulon") + Seurat::NoAxes() + #NoLegend() +
    # theme(plot.title = element_blank())
    ggtitle(gene) + guides(color=guide_colorbar(title="Enrichment"))

  pp <- cowplot::plot_grid(p2, p1, ncol=2)
  return(pp)
}















