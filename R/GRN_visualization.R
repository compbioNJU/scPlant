


#' Heatmaps showing mean regulon activity and TF expression of each cluster.
#'
#' @param SeuratObj Seurat object
#' @param rasMat matrix of regulon activity score(RAS) in each cell
#' @param group.by Name of the metadata column to group cells by (for example, seurat_clusters)
#' @param assay assay to calculate mean expression
#'
#' @importFrom ComplexHeatmap Heatmap draw
#' @importFrom grid unit
#' @importFrom dplyr `%>%`
#'
#' @export
#'
ras_exp_hmp <- function(SeuratObj, rasMat, group.by = "seurat_clusters", assay = 'SCT') {
  # mean regulon activity
  meanRAS <- sapply(split(rownames(SeuratObj@meta.data), SeuratObj@meta.data[[group.by]]),
                    function(cells) colMeans(rasMat[cells, ]))
  meanRAS <- t(scale(t(meanRAS), center = T, scale=T))
  col_fun1 = circlize::colorRamp2(seq(min(meanRAS), max(meanRAS), length.out = 10), c("white", pals::brewer.orrd(9)))

  ph <- pheatmap::pheatmap(meanRAS, cluster_rows = T, cluster_cols = T, silent = T)
  rowOrder <- ph$tree_row$labels[ph$tree_row$order]
  colOrder <- ph$tree_col$labels[ph$tree_col$order]
  # mean TF expression
  Seurat::Idents(SeuratObj) <- SeuratObj@meta.data[[group.by]]
  meanExp <- Seurat::AverageExpression(SeuratObj, assays = assay, slot = 'data')[[assay]]
  meanExp <- meanExp[colnames(rasMat), ] %>% as.matrix()
  meanExp <- t(scale(t(meanExp), center = T, scale=T))
  col_fun2 = circlize::colorRamp2(seq(min(meanExp), max(meanExp), length.out = 10), c("white", pals::brewer.blues(9)))

  ht1 <- Heatmap(meanRAS[rowOrder, colOrder], name = "Regulon Activity", col = col_fun1,
                 column_names_rot = 45, row_names_gp = gpar(fontsize = 1),
                 cluster_rows = F, cluster_columns = F,
                 column_names_gp = gpar(fontsize = 7), column_names_centered = F,  border = TRUE)
  ht2 <- Heatmap(meanExp[rowOrder, colOrder], name = "TF Expression", col = col_fun2,
                 column_names_rot = 45, row_names_gp = gpar(fontsize = 1),
                 cluster_rows = F, cluster_columns = F,
                 column_names_gp = gpar(fontsize = 7), column_names_centered = F,  border = TRUE)
  ht_list = ht1 + ht2
  draw(ht_list, ht_gap = unit(1, "cm"))
}


#' Dimension reduction plot showing regulon activity and TF expression
#'
#' @param SeuratObj Seurat object
#' @param rasMat matrix of regulon activity score(RAS) in each cell
#' @param gene TF to show
#' @param reduction umap or tsne
#'
#' @importFrom Seurat FeaturePlot NoAxes
#' @importFrom ggplot2 ggtitle guides guide_colorbar
#'
#' @export
#'
ras_exp_scatter <- function(SeuratObj, rasMat, gene, reduction = 'umap') {

  p1 <- FeaturePlot(SeuratObj, reduction = reduction, order=TRUE,
                    cols=c('lightgrey', pals::brewer.blues(10)),
                    coord.fixed = F,
                    pt.size=0.2, combine=T, features = gene) + NoAxes() + #NoLegend() +
    # theme(plot.title = element_blank())
  ggtitle(gene) + guides(color=guide_colorbar(title="Expression"))

  SeuratObj$Regulon <- rasMat[rownames(SeuratObj@meta.data), gene]
  p2 <- FeaturePlot(SeuratObj, reduction = reduction, order=TRUE, coord.fixed = F,
                    cols=c("lightgrey", pals::brewer.orrd(10)),
                    pt.size=0.2, combine=T,
                    features="Regulon") + NoAxes() + #NoLegend() +
    # theme(plot.title = element_blank())
  ggtitle(gene) + guides(color=guide_colorbar(title="Regulon Activity"))

  pp <- cowplot::plot_grid(p2, p1, ncol=2)
  return(pp)
}



#' Network diagram showing top regulons of each cluster
#'
#' Network diagram showing top regulons of each cluster according to regulon specificity score(RSS).
#'
#' @param rssMat matrix of regulon specificity score(RSS) in each cluster
#' @param topn number of top regulons to draw
#' @param from.size size of 'cluster' nodes
#' @param to.size size of 'regulon' nodes
#' @param from.label.cex label size of 'cluster' nodes
#' @param to.label.cex label size of 'regulon' nodes
#'
#' @importFrom igraph graph_from_data_frame E V `E<-` `V<-`
#'
#' @export
#'
topRegulons <- function(rssMat, topn = 5,
                        from.size = 7, to.size = 4,
                        from.label.cex = 1, to.label.cex = 0.7) {
  # choose top regulons of each cluster according to regulon specificity score(RSS)
  ll <- lapply(colnames(rssMat), function(x){
    tt <- head(sort(rssMat[, x], decreasing = T), topn)
    data.frame(celltype = x, regulon = names(tt), RSS = tt, stringsAsFactors = F)
  })
  chosen_regulon <- do.call('rbind', ll)
  # network diagram
  edgedf <- chosen_regulon %>% magrittr::set_colnames(c('source', 'target', 'width'))
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



#' Network diagram showing top targets of each regulon
#'
#' Network diagram showing top targets of each regulon according to importance score.
#'
#' @param tf_target TF and targets of each regulon. This information comes from post-processing of pyscenic.
#' @param topn number of top targets to draw
#' @param regulons regulons to show. If NULL, show all regulons.
#' @param tf.size size of TF nodes
#' @param target.size size of target nodes
#' @param tf.label.cex label size of TF nodes
#' @param target.label.cex label size of target nodes
#'
#' @export
#'
toptargets <- function(tf_target, topn = 5, regulons = NULL,
                       tf.size = 7, target.size = 4,
                       tf.label.cex = 1, target.label.cex = 0.7) {
  # choose top targets of each regulon according to importance score
  edges <- tf_target %>% dplyr::filter(TF %in% regulons) %>% dplyr::group_by(TF) %>%
    dplyr::slice_max(order_by = importance_score, n=topn, with_ties = F) %>% dplyr::ungroup()
  edges <- edges[,1:2] %>% magrittr::set_colnames(c("source", "target"))

  net <- graph_from_data_frame(d=edges, directed=T)
  E(net)$width <- 0.2
  E(net)$arrow.size <- .1
  E(net)$arrow.width <- .1
  V(net)$label.cex <- ifelse(V(net)$name %in% unique(edges$source), tf.label.cex, target.label.cex)
  V(net)$label.font <- 4
  V(net)$label.color <- "black"
  V(net)$frame.color <- NA
  V(net)$color <- pals::brewer.greys(9)[3]
  V(net)$size <- ifelse(V(net)$name %in% unique(edges$source), tf.size, target.size)
  V(net)$shape <- 'circle'
  set.seed(123)
  plot(net, layout=igraph::layout_with_kk, vertex.label.degree=0, vertex.label.dist=0)
}



#' Network diagram showing top targets of top regulons of each cluster
#'
#' Each regulator is colored by the cluster that it belongs. Pie means that this regulon is top regulon of multiple clusters.
#'
#' @param rssMat matrix of regulon specificity score(RSS) in each cluster
#' @param tf_target TF and targets of each regulon. This information comes from post-processing of pyscenic.
#' @param Topregulons number of top regulons to draw
#' @param Toptargets number of top targets to draw
#'
#' @export
#'
ToptargetsofTopregulons <- function(rssMat,
                                    tf_target,
                                    Topregulons = 5,
                                    Toptargets = 5) {
  # choose top regulons of each cluster according to regulon specificity score(RSS)
  ll <- lapply(colnames(rssMat), function(x){
    tt <- head(sort(rssMat[, x], decreasing = T), Topregulons)
    data.frame(celltype = x, regulon = names(tt), RSS = tt, stringsAsFactors = F)
  })
  chosen_regulon <- do.call('rbind', ll)
  df <- unique(chosen_regulon[, c('regulon', 'celltype')])
  piedata <- table(df$regulon, df$celltype) %>% as.data.frame()
  piedata <- reshape2::acast(piedata, Var1~Var2, value.var = "Freq")
  piedata <- as.data.frame(piedata[rowSums(piedata)>1, ]) # pie data that we need to draw

  groupColor <- setNames(CellFunTopic::scPalette2(ncol(rssMat)), colnames(rssMat))
  nodecolor <- df %>% dplyr::mutate(celltype = groupColor[celltype]) %>% tibble::deframe()
  piecolor <- groupColor

  edges <- tf_target %>% dplyr::filter(TF %in% unique(chosen_regulon$regulon)) %>%
    dplyr::group_by(TF) %>% dplyr::slice_max(order_by = importance_score, n=Toptargets, with_ties = F) %>% dplyr::ungroup() %>%
    dplyr::select(-importance_score) %>% magrittr::set_colnames(c("source", "target"))

  net <- graph_from_data_frame(d=edges, directed=T)
  E(net)$width <- 0.2
  E(net)$arrow.size <- .1
  E(net)$arrow.width <- .1
  V(net)$label.cex <- ifelse(V(net)$name %in% unique(edges$source), 0.4, 0.2)
  V(net)$label.font <- 4
  V(net)$label.color <- "black"
  V(net)$frame.color <- NA
  V(net)$size <- ifelse(V(net)$name %in% unique(edges$source), 3, 1.5)
  V(net)$shape <- ifelse(V(net)$name %in% rownames(piedata), 'pie', 'circle')
  V(net)$color <- ifelse(V(net)$name %in% unique(edges$source), nodecolor[V(net)$name], pals::brewer.greys(9)[3])
  pie.values <- lapply(V(net)$name, function(x){
    as.numeric(piedata[x, ])
  })
  set.seed(123)
  plot(net, layout=igraph::layout_with_kk, vertex.label.degree=0, vertex.label.dist=0, vertex.pie=pie.values,
       vertex.pie.border=NA, vertex.pie.color=list(piecolor[colnames(piedata)]))
}



#' Dot plot showing top regulons of each cluster according to regulon specificity score(RSS)
#'
#' @param rssMat matrix of regulon specificity score(RSS) in each cluster
#' @param cluster which cluster to show
#' @param topn number of top regulons to draw
#'
#' @importFrom ggplot2 geom_point ylab xlab theme_bw theme element_blank element_line element_text
#'
#' @export
#'
SpecificityRank <- function(rssMat, cluster, topn = 10) {
  data <- tibble::enframe(sort(rssMat[, as.character(cluster)], decreasing = T)) %>%
    magrittr::set_colnames(c('regulon', 'RSS')) %>%
    tibble::rowid_to_column("index")

  data$pt.col <- ifelse(data$index <= topn, "#007D9B", "#BECEE3")
  data <- head(data, n=200)
  data.label <- head(data, n=topn)

  pp <- ggplot(data, aes(index, RSS)) +
    geom_point(size=3, color=data$pt.col) +
    ggrepel::geom_text_repel(inherit.aes = FALSE, data = data.label, aes(index, RSS, label=regulon), size=4,
                             max.overlaps = 30) +
    ggtitle(as.character(cluster)) + ylab("Specificity score") + xlab('Regulons') +
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



