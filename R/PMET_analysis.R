
#' pre-process original PMET output and filter with PPI information
#'
#' @param PMETresult path of original PMET output
#' @param PPIinfo protein-protein interaction (PPI) information
#' @param species species; currently \code{Arabidopsis thaliana} or \code{Oryza sativa} or \code{Zea mays} are supported.
#'
#' @importFrom magrittr `%>%` `%<>%`
#' @importFrom utils read.table
#'
#' @export
#'
#' @examples
#' \dontrun{
#' PMETresult <- processPMET(PMETresult = "output.txt")
#' }
#'
processPMET <- function(PMETresult, PPIinfo = NULL, species = 'Arabidopsis thaliana') {
  if (species == 'Arabidopsis thaliana') {
    motif2TF <- motif2TF_At
  } else if (species == 'Oryza sativa') {
    motif2TF <- motif2TF_Os
  } else if (species == 'Zea mays') {
    motif2TF <- motif2TF_Zm
  }
  # motif-to-TF mapping
  PMETresult <- read.table(file = PMETresult, header = T, sep = "\t", stringsAsFactors = F)
  PMETresult %<>% dplyr::mutate(TF1 = unname(motif2TF[Motif1]), TF2 = unname(motif2TF[Motif2]))
  if (!is.null(PPIinfo)) {
    PPIinfo <- PPIinfo[, 2:3]
    # only retain TF-TF pairs with protein-protein interaction (PPI) evidence
    PPIinfo <- data.table::rbindlist(list(PPIinfo, PPIinfo[, 2:1]), use.names = F) %>% unique()
    ii <- paste(PMETresult$TF1, PMETresult$TF2, sep = "_") %in% paste(PPIinfo$INTERACTOR_A, PPIinfo$INTERACTOR_B, sep = "_")
    PMETresult <- PMETresult[ii,]
    # PMETresult <- dplyr::semi_join(PMETresult, PPIinfo, by = c("TF1" = "INTERACTOR_A", "TF2" = "INTERACTOR_B")) # code of same effect
  }
  PMETresult %<>% dplyr::filter(TF1 != TF2)
  # Considering the one-to-many mapping between TF and motif, we only retain TF pair with the highest enrichment out of replicated pairs.
  PMETresult$pair <- apply(PMETresult, 1, FUN = function(x){paste(sort(c(x[5], x[6])), collapse = "_")})
  PMETresult %<>% dplyr::group_by(Module, pair) %>% dplyr::slice_min(order_by = AdjustedpvalueBH, n=1, with_ties=F) %>%
    dplyr::ungroup()
  PMETresult %<>% dplyr::mutate(logFDR=log1p(-log10(AdjustedpvalueBH)))
  return(PMETresult)
}


#' Heatmap of PMET result
#'
#' @param PMETresult PMET result
#' @param topn show top highly enriched pairs of each cluster
#' @param rownames.size font size of row names
#' @param colnames.size font size of column names
#'
#' @importFrom ComplexHeatmap Heatmap draw
#'
#' @export
#'
PMEThmp <- function(PMETresult, topn = NULL, rownames.size = 5, colnames.size = 5) {

  mm <- reshape2::acast(PMETresult, Module~pair, value.var='logFDR')
  mm <- mm[, colSums(mm) > 0]
  if (!is.null(topn)) {
    pairs <- PMETresult %>% dplyr::group_by(Module) %>% dplyr::slice_max(order_by = logFDR, n = topn, with_ties = F) %>%
      dplyr::ungroup() %>% dplyr::pull(pair) %>% unique
    mm <- mm[, pairs]
  }
  mm <- t(mm)
  col_fun = circlize::colorRamp2(seq(min(mm), max(mm), length.out = 10), c("white", pals::brewer.greens(9)))
  mm <- t(scale(t(mm), center = T, scale=T))
  ht <- Heatmap(mm, name = "log1p(-log10(P.adjust))", col = col_fun,
                column_names_rot = 45, row_names_gp = gpar(fontsize = rownames.size),
                column_names_gp = gpar(fontsize = colnames.size),
                cluster_rows = T, cluster_columns = T, column_names_centered = T, border = TRUE)
  draw(ht)
}


#' Draw triangle heatmap of PMET result
#'
#' @param PMETresult PMET result
#' @param clus cluster
#'
#' @importFrom grDevices colorRampPalette
#' @importFrom ggplot2 geom_tile scale_fill_gradientn theme_classic coord_fixed scale_x_discrete labs
#'
#' @export
#'
triHmp <- function(PMETresult, clus) {
  df <- PMETresult[PMETresult$Module == clus, ]
  # only show top pairs of high enrichment
  dff <- df %>% dplyr::slice_max(order_by = logFDR, prop=0.1)
  tfs <- unique(c(dff$TF1, dff$TF2))
  dff <- df %>% dplyr::filter((TF1 %in% tfs) & (TF2 %in% tfs))
  dff <- data.table::rbindlist(list(dff[, c("Module",  "TF1","TF2", "logFDR")],
                                    dff[, c("Module",  "TF2","TF1", "logFDR")]), use.names = F) %>% unique()
  mat <- reshape2::acast(dff, TF1~TF2, value.var = "logFDR")
  mat[is.na(mat)] <- 0

  colss <- c("#FFFFFF", colorRampPalette(c("#FFFFCC", "#EF3B2C", "#9932CC", "#000000"))(100))
  pp <- pheatmap::pheatmap(mat, color = colss, silent=T)
  mat <-mat[pp$tree_row$order, pp$tree_row$order]
  mat[lower.tri(mat)] <- NA
  melted_mat <- reshape2::melt(mat, na.rm = TRUE)

  pp <- ggplot(data = melted_mat, aes(Var2, Var1, fill = value))+
    geom_tile() + coord_fixed(ratio=1) +
    scale_fill_gradientn(colours = colss, name="log1p(-log10(padj))") +
    theme_classic()+
    theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 5, hjust = 1),
          axis.text.y = element_text(size = 5),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          axis.line = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    scale_x_discrete(limits=rev(levels(melted_mat$Var2))) +
    labs(title = as.character(clus))
  return(pp)
}


#' Network diagram showing top pairs of a cluster.
#'
#' @param PMETresult PMET result
#' @param clus cluster
#' @param topn number of top highly enriched pairs of the cluster to draw, 20 as default.
#' @param nodeColor color of nodes
#'
#' @importFrom igraph graph_from_data_frame incident_edges
#'
#' @export
#'
topPairsNet <- function(PMETresult, clus, topn = 20, nodeColor = "#81BADA") {
  edges <- PMETresult %>% dplyr::filter(Module == clus) %>%
    dplyr::slice_max(order_by = logFDR, n=topn, with_ties=F) %>%
    dplyr::select(TF1, TF2, logFDR)
  colnames(edges) <- c('source', 'target', 'width')
  edges$width <- scales::rescale(edges$width, to = c(1, 4))
  net <- graph_from_data_frame(edges, directed = F)
  V(net)$color <- nodeColor
  V(net)$label.cex <- 0.7
  V(net)$label.font <- 4
  V(net)$label.color <- "black"
  inc.edges <- incident_edges(net, V(net), mode="all")
  V(net)$size <- scales::rescale(unlist(lapply(inc.edges, length)), to = c(9, 15))
  set.seed(123)
  plot(net, vertex.label.dist=2)
}





