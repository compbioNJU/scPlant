#' Run monocle2 on Seurat object
#'
#' @param SeuratObj A Seurat object
#' @param max_components the dimensionality of the reduced space. 2 as default.
#' @param root_clusters clusters selected as root
#' @param point_size point size of scatter plot
#' @param BEAM_analysis perform BEAM analysis or not.
#' @param BEAM_branch_point The ID of the branch point to analyze.
#' @param save_plot save plots to pdf or not.
#' @param return_object return CellDataSet object or not.
#'
#' @importFrom grDevices dev.off pdf
#'
#' @return if parameter \code{return_object = TRUE}, CellDataSet object and differential genes across pseudotime will be returned.
#'
#' @export
#'
RunMonocle2 <- function(SeuratObj,
                        max_components = 2,
                        root_clusters = NULL,
                        point_size = 0.5,
                        BEAM_analysis = FALSE,
                        BEAM_branch_point = 1,
                        save_plot = TRUE,
                        return_object = TRUE) {

  #### create CellDataSet object
  data <- Seurat::GetAssayData(SeuratObj, slot = "counts")
  pd <- SeuratObj@meta.data
  pd$cds_cluster <- unname(Seurat::Idents(SeuratObj)[rownames(pd)])
  fData <- data.frame(gene_short_name = rownames(data), geneID=rownames(data), row.names = rownames(data))
  mycds <- monocle::newCellDataSet(data,
                                   phenoData = new('AnnotatedDataFrame', data = pd),
                                   featureData = new('AnnotatedDataFrame', data = fData),
                                   expressionFamily = negbinomial.size())
  #### run Monocle2 pipeline
  mycds <- estimateSizeFactors(mycds)
  mycds <- estimateDispersions(mycds, relative_expr = TRUE)
  # we use marker genes obtained from Seurat pipeline
  diff.genes <- slot(object = SeuratObj, name = 'misc')[["Allmarkers"]]
  sig_diff.genes <- subset(diff.genes, p_val_adj<0.01 & abs(avg_log2FC)>0.5)$gene
  sig_diff.genes <- unique(as.character(sig_diff.genes))

  mycds <- monocle::setOrderingFilter(mycds, sig_diff.genes)
  mycds <- monocle::reduceDimension(mycds, max_components = max_components, method = 'DDRTree')
  mycds <- monocle::orderCells(mycds)
  # select root state
  if (!is.null(root_clusters)) {
    root_state <- function(mycds){
      if (length(unique(pData(mycds)$State)) > 1){
        R_counts <- table(pData(mycds)$State, pData(mycds)$cds_cluster)[,as.character(root_clusters), drop=F] %>% rowSums()
        # return(as.numeric(names(R_counts)[which(R_counts == max(R_counts))]))
        return(as.numeric(names(R_counts)[which.max(R_counts)]))
      } else {
        return(1)
      }
    }
    mycds <- orderCells(mycds, root_state = root_state(mycds))
  }

  #### differential genes across pseudotime
  diff_test <- monocle::differentialGeneTest(mycds[sig_diff.genes,], cores = 1,
                                             fullModelFormulaStr = "~sm.ns(Pseudotime)")

  #### BEAM analysis
  if (BEAM_analysis) {
    disp_table <- monocle::dispersionTable(mycds)
    disp.genes <- subset(disp_table, mean_expression >= 0.5 & dispersion_empirical >= 1*dispersion_fit)
    disp.genes <- as.character(disp.genes$gene_id)
    mycds_sub <- mycds[disp.genes,]
    beam_result <- monocle::BEAM(mycds_sub, branch_point = BEAM_branch_point, cores = 1)
  }

  #### save plots
  if (save_plot) {
    if (!dir.exists(paths = "./output/RunMonocle2")) {
      dir.create("./output/RunMonocle2", recursive = TRUE)
    }
    # visualization
    p1 <- monocle::plot_cell_trajectory(mycds, color_by = "State", cell_size = point_size)
    p2 <- monocle::plot_cell_trajectory(mycds, color_by = "cds_cluster", cell_size = point_size)
    p3 <- monocle::plot_cell_trajectory(mycds, color_by = "Pseudotime", cell_size = point_size)
    p4 <- monocle::plot_cell_trajectory(mycds, color_by = "State", cell_size = point_size) + facet_wrap(~State, ncol = 2)
    p5 <- monocle::plot_cell_trajectory(mycds, color_by = "cds_cluster", cell_size = point_size) + facet_wrap(~cds_cluster, ncol = 4)
    sig_gene_names <- rownames(subset(diff_test, qval < 0.01))
    p6 <- monocle::plot_pseudotime_heatmap(mycds[sig_gene_names,], num_clusters=3, use_gene_short_name = TRUE,
                                           show_rownames=T, return_heatmap=T)
    pdf(file = "./output/RunMonocle2/cell_trajectory.pdf", width = 15, height = 13)
    print(p1)
    print(p2)
    print(p3)
    dev.off()
    pdf(file = "./output/RunMonocle2/cell_trajectory_facet.pdf", width = 10, height = 20)
    print(p4)
    print(p5)
    dev.off()
    pdf(file = "./output/RunMonocle2/pseudotime_heatmap.pdf", width = 10, height = 17)
    print(p6)
    dev.off()
    if (BEAM_analysis) {
      beam_res <- beam_result[order(beam_result$qval),]
      beam_res <- beam_res[,c("gene_short_name", "pval", "qval")]
      mycds_sub_beam <- mycds_sub[rownames(subset(beam_res, qval < 1e-4)),]
      p7 <- monocle::plot_genes_branched_heatmap(mycds_sub_beam,  branch_point = BEAM_branch_point, num_clusters = 3,
                                                 show_rownames = T, use_gene_short_name = TRUE, return_heatmap = TRUE)
      pdf(file = "./output/RunMonocle2/branched_heatmap.pdf", width = 10, height = 17)
      print(p7)
      dev.off()
    }
  }
  if (return_object) {
    if (BEAM_analysis) {
      return(list(CellDataSet=mycds, Gene_acrossPseudotime=diff_test, beam_result=beam_result))
    } else {
      return(list(CellDataSet=mycds, Gene_acrossPseudotime=diff_test))
    }
  }
}
















