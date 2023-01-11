



################### 调控网络部分的画图函数
topRegulons <- function(rssMat, topn = 5) {
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
  V(net)$label.cex <- ifelse(V(net)$name %in% celltypes, 1, 0.7)
  V(net)$label.font <- ifelse(V(net)$name %in% celltypes, 2, 3)
  V(net)$label.color <- "black"
  V(net)$frame.color <- NA
  V(net)$size <- ifelse(V(net)$name %in% celltypes, 7, 3)
  V(net)$shape <- "circle"
  V(net)$color <- ifelse(V(net)$name %in% celltypes, nodeColor[V(net)$name], pals::brewer.greys(9)[3])
  V(net)$label.dist <- ifelse(V(net)$name %in% celltypes, 0, 0.5)
  set.seed(123)
  plot(net, vertex.label.degree=0, vertex.label.dist=0)
}


toptargets <- function(tf_target, topn = 5, regulons = NULL) {
  # choose top targets of each regulon according to importance score
  edges <- tf_target %>% dplyr::filter(TF %in% regulons) %>% dplyr::group_by(TF) %>% 
    dplyr::slice_max(order_by = importance_score, n=topn, with_ties = F) %>% dplyr::ungroup()
  edges <- edges[,1:2] %>% magrittr::set_colnames(c("source", "target"))
  
  net <- graph_from_data_frame(d=edges, directed=T) 
  E(net)$width <- 0.2
  E(net)$arrow.size <- .1
  E(net)$arrow.width <- .1
  V(net)$label.cex <- ifelse(V(net)$name %in% unique(edges$source), 1, 0.7)
  V(net)$label.font <- ifelse(V(net)$name %in% unique(edges$source), 2, 3)
  V(net)$label.color <- "black"
  V(net)$frame.color <- NA
  V(net)$color <- pals::brewer.greys(9)[3]
  V(net)$size <- ifelse(V(net)$name %in% unique(edges$source), 8, 5)
  V(net)$shape <- 'circle'
  set.seed(123)
  plot(net, layout=layout_with_kk, vertex.label.degree=0, vertex.label.dist=0)
}



ras_exp_hmp2 <- function(meanRAS, meanExp, regulons) {
  meanRAS <- meanRAS[regulons, ]
  # mean TF expression
  meanExp <- meanExp[regulons, ] %>% as.matrix()
  meanExp <- t(scale(t(meanExp), center = T, scale=T))
  if (length(regulons) > 1) {
    ph <- pheatmap::pheatmap(meanRAS, cluster_rows = T, cluster_cols = F, silent = T)
    rowOrder <- ph$tree_row$labels[ph$tree_row$order]
    meanRAS <- meanRAS[rowOrder, ]
    meanExp <- meanExp[rowOrder, ]
  }
 
  col_fun1 = circlize::colorRamp2(seq(min(meanRAS), max(meanRAS), length.out = 10), c("white", pals::brewer.orrd(9)))
  col_fun2 = circlize::colorRamp2(seq(min(meanExp), max(meanExp), length.out = 10), c("white", pals::brewer.blues(9)))
  
  ht1 <- ComplexHeatmap::Heatmap(meanRAS, name = "Regulon Activity", col = col_fun1, 
                                 column_names_rot = 45, row_names_gp = gpar(fontsize = 10),
                                 cluster_rows = F, cluster_columns = F, 
                                 column_names_gp = gpar(fontsize = 12), column_names_centered = T,  border = TRUE)
  ht2 <- ComplexHeatmap::Heatmap(meanExp, name = "TF Expression", col = col_fun2, 
                                 column_names_rot = 45, row_names_gp = gpar(fontsize = 10),
                                 cluster_rows = F, cluster_columns = F, 
                                 column_names_gp = gpar(fontsize = 12), column_names_centered = T,  border = TRUE)
  ht_list = ht1 + ht2
  ComplexHeatmap::draw(ht_list, ht_gap = unit(1, "cm"))
}


########### 发现用pheatmap画的topic-cluster热图经常显示不出来，改用ComplexHeatmap画

cluster_topic_hmp2 <- function(ldaOut) 
{
  mm <- posterior(ldaOut)$topics
  colnames(mm) <- paste0("Topic ", colnames(mm))
  col_fun = circlize::colorRamp2(seq(min(mm), max(mm), length.out = 10), c("white", pals::brewer.orrd(9)))
  ht <- ComplexHeatmap::Heatmap(mm, name = "Topic probability", col = col_fun, 
                                 column_names_rot = 45, row_names_gp = gpar(fontsize = 10),
                                 cluster_rows = T, cluster_columns = T, rect_gp = gpar(col = "grey"),
                                 column_names_gp = gpar(fontsize = 9), column_names_centered = T,  border = TRUE)
  ComplexHeatmap::draw(ht)
}















