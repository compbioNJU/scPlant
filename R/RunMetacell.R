

#' Run metacell on Seurat object
#'
#' The MetaCell R package facilitates analysis of single cell RNA-seq UMI matrices by computing partitions
#' of a cell similarity graph into small (~20-200 typically) homogeneous groups of cells which are defined as
#' metacells (MCs). The derived MCs are then used for building different representations of the data,
#' allowing matrix or 2D graph visualization forming a basis for analysis of cell types, subtypes,
#' transcriptional gradients, cell-cycle variation, gene modules and their regulatory models and more.
#' More details on the usage of the MetaCell pipeline is available in the package vignettes \url{https://tanaylab.github.io/metacell/index.html}.
#' Note that metacell package is tested on linux and macbooks, and is currently not compatible on Windows.
#'
#' @param SeuratObj Seurat object
#' @param species species
#'
#' @export
#'
RunMetacell <- function(SeuratObj, species) {

  #  initialize a database, linking the package to directory that stores all your objects
  if(!dir.exists("metacelldb")) dir.create("metacelldb/")
  metacell::scdb_init("metacelldb/", force_reinit=T)
  if ("SCT" %in% names(SeuratObj@assays)) {
    SeuratObj <- Seurat::as.SingleCellExperiment(SeuratObj, assay = "SCT")
  } else {
    SeuratObj <- Seurat::as.SingleCellExperiment(SeuratObj)
  }
  mat = metacell::scm_import_sce_to_mat(SeuratObj)
  # link the package to a figure directory
  if(!dir.exists("metacellfigs")) dir.create("metacellfigs/")
  metacell::scfigs_init("metacellfigs/")
  # filtering the UMI matrix
  metacell::mcell_plot_umis_per_cell("test")
  # mitochondrial genes and immunoglobulin genes
  MTgenes <- getMTgenes(species)
  nms = c(rownames(mat@mat), rownames(mat@ignore_gmat))
  ig_genes = c(grep("^IGJ", nms, v=T),
               grep("^IGH",nms,v=T),
               grep("^IGK", nms, v=T),
               grep("^IGL", nms, v=T))
  bad_genes = unique(c(MTgenes,"NEAT1","TMSB4X", "TMSB10", ig_genes))
  metacell::mcell_mat_ignore_genes(new_mat_id="test", mat_id="test", bad_genes, reverse=F)
  metacell::mcell_mat_ignore_small_cells("test", "test", 800)

  # Selecting feature genes
  metacell::mcell_add_gene_stat(gstat_id="test", mat_id="test", force=T)
  metacell::mcell_gset_filter_varmean(gset_id="test_feats", gstat_id="test", T_vm=0.08, force_new=T)
  metacell::mcell_gset_filter_cov(gset_id = "test_feats", gstat_id="test", T_tot=100, T_top3=2)
  metacell::mcell_plot_gstats(gstat_id="test", gset_id="test_feats")
  # Building the balanced cell graph
  metacell::mcell_add_cgraph_from_mat_bknn(mat_id="test",
                                 gset_id = "test_feats",
                                 graph_id="test_graph",
                                 K=100,
                                 dsamp=T)
  # Resampling and generating the co-clustering graph
  metacell::mcell_coclust_from_graph_resamp(coc_id="test_coc500",
                                  graph_id="test_graph",
                                  min_mc_size=20,
                                  p_resamp=0.75, n_resamp=500)
  metacell::mcell_mc_from_coclust_balanced(coc_id="test_coc500",
                                 mat_id= "test",
                                 mc_id= "test_mc",
                                 K=30, min_mc_size=30, alpha=2)
  # Removing outlier cells
  metacell::mcell_plot_outlier_heatmap(mc_id="test_mc", mat_id = "test", T_lfc=3)
  metacell::mcell_mc_split_filt(new_mc_id="test_mc_f",
                      mc_id="test_mc",
                      mat_id="test",
                      T_lfc=3, plot_mats=F)
  # Selecting markers and coloring metacells
  metacell::mcell_gset_from_mc_markers(gset_id="test_markers", mc_id="test_mc_f")
  metacell::mc_colorize_default("test_mc_f")
  # Creating a heatmap of genes and metacells
  metacell::mcell_mc_plot_marks(mc_id="test_mc_f", gset_id="test_markers", mat_id="test")
  # Projecting metacells and cells in 2D
  metacell::mcell_mc2d_force_knn(mc2d_id="test_2dproj",mc_id="test_mc_f", graph_id="test_graph")
  tgconfig::set_param("mcell_mc2d_height",1000, "metacell")
  tgconfig::set_param("mcell_mc2d_width",1000, "metacell")
  metacell::mcell_mc2d_plot(mc2d_id="test_2dproj")
  # Visualizing the MC confusion matrix
  mc_hc = metacell::mcell_mc_hclust_confu(mc_id="test_mc_f",
                                graph_id="test_graph")
  mc_sup = metacell::mcell_mc_hierarchy(mc_id="test_mc_f",
                              mc_hc=mc_hc, T_gap=0.04)
  metacell::mcell_mc_plot_hierarchy(mc_id="test_mc_f",
                          graph_id="test_graph",
                          mc_order=mc_hc$order,
                          sup_mc = mc_sup,
                          width=2800, heigh=2000, min_nmc=2)


}










