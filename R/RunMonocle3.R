

#' Run monocle3 on Seurat object
#'
#' @param SeuratObj A Seurat object
#' @param batch one column of \code{SeuratObj@meta.data}, indicating batch effect.
#' @param reduction_method reduction method, UMAP as default.
#' @param keep_orig_reduction keep original reduction dimension of Seurat object.
#' @param root_clusters seurat clusters used to help identify the root.
#' @param save_plot save plots to pdf or not.
#' @param return_object return CellDataSet object or not.
#' 
#' @importFrom grDevices dev.off pdf
#'
#' @return if \code{return_object = TRUE}, CellDataSet object, differential genes across pseudotime, 
#' differential genes between clusters, gene modules will be returned.
#' @export
#'
RunMonocle3 <- function(SeuratObj, 
                        batch = NULL,
                        reduction_method = "UMAP",
                        keep_orig_reduction = TRUE,
                        root_clusters = NULL, 
                        save_plot = TRUE,
                        return_object = TRUE) {
  
  #### create CellDataSet object
  data <- Seurat::GetAssayData(SeuratObj, assay = 'RNA', slot = 'counts')
  cell_metadata <- SeuratObj@meta.data
  gene_annotation <- data.frame(gene_short_name = rownames(data))
  rownames(gene_annotation) <- rownames(data)
  cds <- monocle3::new_cell_data_set(data,
                                     cell_metadata = cell_metadata,
                                     gene_metadata = gene_annotation)
  
  #### Pre-process the data
  cds <- monocle3::preprocess_cds(cds, num_dim = 50)
  
  #### remove batch effects
  if (!is.null(batch)) {
    cds <- monocle3::align_cds(cds, alignment_group = batch)
  }
  
  #### Reduce dimensionality
  cds <- monocle3::reduce_dimension(cds, reduction_method=reduction_method)
  if (keep_orig_reduction) {
    cds.embed <- cds@int_colData$reducedDims[[reduction_method]]
    seurat.embed <- Seurat::Embeddings(SeuratObj, reduction=tolower(reduction_method))
    cds@int_colData$reducedDims[[reduction_method]] <- seurat.embed[rownames(cds.embed), ]
  }
  
  #### Group cells into clusters and partitions
  cds <- monocle3::cluster_cells(cds, resolution=1e-5)
  
  #### Learn the trajectory graph
  cds <- monocle3::learn_graph(cds)
  
  # select root
  if (!is.null(root_clusters)) {
    # a helper function to identify the root principal points
    get_earliest_principal_node <- function(cds, clus = root_clusters){
      cell_ids <- which(colData(cds)[, "seurat_clusters"] %in% clus)
      
      closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
      closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
      root_pr_nodes <-
        igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
      
      root_pr_nodes
    }
    cds <- monocle3::order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))
  }
  
  #### Differential expression analysis
  #### Graph-autocorrelation analysis for comparing clusters
  deg_clusters <- monocle3::graph_test(cds, neighbor_graph="knn", cores=1)
  # Finding modules of co-regulated genes
  pr_deg_ids <- row.names(subset(deg_clusters, q_value < 0.05))
  deg_clusters_module <- monocle3::find_gene_modules(cds[pr_deg_ids,], resolution=1e-2)
  if (save_plot) {
    cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
                                    cell_group=partitions(cds)[colnames(cds)])
    agg_mat <- monocle3::aggregate_gene_expression(cds, deg_clusters_module, cell_group_df)
    row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
    colnames(agg_mat) <- stringr::str_c("Partition ", colnames(agg_mat))
    
    pheatmap::pheatmap(agg_mat, cluster_rows=TRUE, cluster_cols=TRUE,
                       scale="column", clustering_method="ward.D2",
                       fontsize=6, filename = './output/plots/RunMonocle3/deg_clusters_module.pdf')
  }
  
  #### Finding genes that change as a function of pseudotime
  deg_pseudotime <- monocle3::graph_test(cds, neighbor_graph="principal_graph", cores=1)
  # collect the trajectory-variable genes into modules:
  pr_deg_ids <- row.names(subset(deg_pseudotime, q_value < 0.01))
  deg_pseudotime_module <- monocle3::find_gene_modules(cds[pr_deg_ids,], resolution=c(10^seq(-6,-1)))
  if (save_plot) {
    cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
                                    cell_group=colData(cds)$seurat_clusters)
    agg_mat <- monocle3::aggregate_gene_expression(cds, deg_pseudotime_module, cell_group_df)
    row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
    
    pheatmap::pheatmap(agg_mat, scale="column", clustering_method="ward.D2", 
                       filename = './output/plots/RunMonocle3/deg_pseudotime_module.pdf')
  }
  
  #### Find marker genes expressed by each cluster
  marker_test_res <- monocle3::top_markers(cds, group_cells_by="partition", 
                                 reference_cells=1000, cores=1)
  if (save_plot) {
    top_specific_markers <- marker_test_res %>%
      dplyr::filter(fraction_expressing >= 0.10) %>%
      dplyr::group_by(cell_group) %>%
      dplyr::top_n(3, pseudo_R2)
    top_specific_marker_ids <- unique(top_specific_markers %>% dplyr::pull(gene_id))
    pdf('./output/plots/RunMonocle3/plot_genes_by_group.pdf')
    monocle3::plot_genes_by_group(cds,
                                  top_specific_marker_ids,
                                  group_cells_by="partition",
                                  ordering_type="cluster_row_col",
                                  max.size=3)
    dev.off()
  }
  
  
  #### save plots
  if (save_plot) {
    if (!dir.exists(paths = "./output/plots/RunMonocle3")) {
      dir.create("./output/plots/RunMonocle3", recursive = TRUE)
    }
    monocle3::save_monocle_objects(cds=cds, directory_path='./output/plots/RunMonocle3/cds_objects')
    # visualization
    pdf('./output/plots/RunMonocle3/graph.pdf')
    monocle3::plot_cells(cds,
                         color_cells_by = "seurat_clusters",
                         label_groups_by_cluster=FALSE,
                         label_leaves=T,
                         label_branch_points=T)
    dev.off()
    
    pdf('./output/plots/RunMonocle3/pseudotime.pdf', width = 10, height = 8)
    monocle3::plot_cells(cds,
                         color_cells_by = "pseudotime",
                         label_cell_groups=F,
                         label_leaves=T,
                         label_branch_points=T,
                         graph_label_size=1.5)
    dev.off()
    
    pdf('./output/plots/RunMonocle3/genes_expression.pdf')
    monocle3::plot_cells(cds, genes=head(rownames(deg_pseudotime)[order(deg_pseudotime$q_value)], 4),
                         show_trajectory_graph=FALSE,
                         label_cell_groups=FALSE,
                         label_leaves=FALSE)
    dev.off()
    
    # give a clearer view of a gene's dynamics along a single path
    genes <- head(rownames(deg_pseudotime)[order(deg_pseudotime$morans_I, decreasing = T)], 3)
    sub_cds <- cds[rowData(cds)$gene_short_name %in% genes, ]
    pdf('./output/plots/RunMonocle3/plot_genes_in_pseudotime.pdf')
    monocle3::plot_genes_in_pseudotime(sub_cds,
                                       color_cells_by="seurat_clusters",
                                       min_expr=0.5)
    dev.off()
    
    
    pdf('./output/plots/RunMonocle3/plots.pdf')
    monocle3::plot_cells(cds)
    monocle3::plot_cells(cds, color_cells_by="seurat_clusters")
    monocle3::plot_cells(cds, color_cells_by="partition", group_cells_by="partition")
    dev.off()
    
  }
  
  if (return_object) {
    return(list(CellDataSet=cds, 
                deg_clusters=deg_clusters, 
                deg_clusters_module=deg_clusters_module,
                deg_pseudotime=deg_pseudotime,
                deg_pseudotime_module=deg_pseudotime_module,
                marker_test_res=marker_test_res))
  }
  
}





















