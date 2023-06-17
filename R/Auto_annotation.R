

#' Automatic annotation of Seurat object using SingleR package.
#'
#' Predict the cell type of Seurat object using SingleR package.
#'
#' @param SeuratObj Seurat object to annotate
#' @param ref Reference dataset, bulk RNA-seq or microarray data or single-cell data,
#' \code{SummarizedExperiment} or \code{SingleCellExperiment} object containing a matrix of log-expression values
#' with labels stored as \code{ref$label}.
#' @param ref_type type of reference dataset, \code{"bulk"} (bulk RNA-seq) or \code{"single-cell"} (single-cell data) are supported.
#'
#' @importFrom Seurat GetAssayData
#'
#' @return Seurat object. The predicted cell type label is stored in \code{SeuratObj$predicted_label}.
#' @export
#'
AutoAnnotate_SingleR <- function(SeuratObj, ref, ref_type = c("bulk", "single-cell")) {

  ref_type <- match.arg(ref_type)
  norm_count <- GetAssayData(SeuratObj, slot="data")
  if (ref_type == "bulk") {
    pred <- SingleR::SingleR(test = norm_count,
                             ref = ref,
                             labels = ref$label)
  } else {
    pred <- SingleR::SingleR(test = norm_count,
                             ref = ref,
                             labels = ref$label,
                             de.method="wilcox")
  }
  SeuratObj$predicted_label <- pred$labels
  return(SeuratObj)
}



#' Automatic annotation of Seurat object using CellAssign package.
#'
#' @param SeuratObj Seurat object to annotate
#' @param markerDF A data frame with two columns (celltype, gene), containing the markers of each cell type.
#'
#' @return Seurat object. The predicted cell type label is stored in \code{SeuratObj$predicted_label}.
#' @export
#'
AutoAnnotate_CellAssign <- function(SeuratObj, markerDF) {
  marker_list <- split(markerDF$gene, markerDF$celltype)
  # turn marker list into the binary marker by cell type matrix
  marker_mat <- cellassign::marker_list_to_mat(marker_list)

  sceset <- Seurat::as.SingleCellExperiment(SeuratObj, assay = Seurat::DefaultAssay(SeuratObj))
  qclust <- scran::quickCluster(sceset, min.size = 30)
  sceset <- scran::computeSumFactors(sceset, sizes = 15, clusters = qclust)

  oo <- intersect(rownames(marker_mat), rownames(sceset))
  marker_mat <- marker_mat[oo, ]
  sceset <- sceset[oo,]
  s <- SingleCellExperiment::sizeFactors(sceset)
  fit <- cellassign::cellassign(exprs_obj = sceset,
                                marker_gene_info = marker_mat,
                                s = s,
                                learning_rate = 1e-2,
                                shrinkage = TRUE,
                                verbose = FALSE)
  SeuratObj$predicted_label <- cellassign::celltypes(fit)
  return(SeuratObj)
}


#' Automatic annotation of Seurat object using Celaref package.
#'
#' @param SeuratObj Seurat object to annotate
#' @param counts_ref  counts matrix of reference data
#' @param cellinfo_ref data frame of cell information of reference data, the first column as cell ID,
#' the second column as cluster ("Cluster" as column name)
#' @param plot produce plots or not
#'
#' @return data frame of prediction result
#' @export
#'
AutoAnnotate_Celaref <- function(SeuratObj, counts_ref, cellinfo_ref, plot = TRUE) {
  # Load data
  counts_query <- GetAssayData(SeuratObj, slot="data")
  cellinfo_query <- Seurat::Idents(SeuratObj) %>% tibble::enframe() %>% magrittr::set_colnames(c("CellId", "Cluster"))
  query_se <- celaref::load_se_from_tables(counts_matrix   = counts_query,
                                  cell_info_table = cellinfo_query,
                                  group_col_name  = "Cluster")
  ref_se <- celaref::load_se_from_tables(counts_matrix   = counts_ref,
                                cell_info_table = cellinfo_ref,
                                group_col_name  = "Cluster")
  # Filter data
  query_se <- celaref::trim_small_groups_and_low_expression_genes(query_se)
  ref_se <- celaref::trim_small_groups_and_low_expression_genes(ref_se)
  # Setup within-experiment differential expression
  de_table.ref   <- celaref::contrast_each_group_to_the_rest(ref_se, dataset_name="ref")
  de_table.query <- celaref::contrast_each_group_to_the_rest(query_se, dataset_name="query")
  # Plot
  if (plot) {
    pdf(file = "AutoAnnotate_Celaref.pdf")
    print(celaref::make_ranking_violin_plot(de_table.test=de_table.query, de_table.ref=de_table.ref))
    dev.off()
  }
  # And get group labels
  result <- celaref::make_ref_similarity_names(de_table.query, de_table.ref)
  return(result)
}



#' Automatic annotation of Seurat object using Garnett package.
#'
#' @param SeuratObj Seurat object to annotate
#' @param marker_file_path see https://cole-trapnell-lab.github.io/garnett/docs/ to construct a marker file
#'
#' @importFrom Seurat GetAssayData Idents
#'
#' @return Seurat object. The predicted cell type label is stored in \code{SeuratObj$predicted_label}.
#' @export
#'
AutoAnnotate_Garnett <- function(SeuratObj,
                                 marker_file_path) {
  # create CellDataSet object
  data <- Seurat::GetAssayData(SeuratObj, slot = "counts")
  pd <- SeuratObj@meta.data
  pd$garnett_cluster <- unname(Seurat::Idents(SeuratObj)[rownames(pd)])
  fData <- data.frame(gene_short_name = rownames(data), geneID=rownames(data), row.names = rownames(data))
  mycds <- monocle::newCellDataSet(as(data, "dgCMatrix"),
                          phenoData = new('AnnotatedDataFrame', data = pd),
                          featureData = new('AnnotatedDataFrame', data = fData),
                          expressionFamily = negbinomial.size())
  mycds <- estimateSizeFactors(mycds)
  # Train the classifier
  set.seed(260)
  classifier <- garnett::train_cell_classifier(cds = mycds,
                                               marker_file = marker_file_path,
                                               db = "none",
                                               num_unknown = 5)
  # Classifying your cells
  mycds <- garnett::classify_cells(mycds, classifier,
                                   db = "none",
                                   cluster_extend = F)
  SeuratObj$predicted_label <- pData(mycds)[Seurat::Cells(SeuratObj), 'cell_type']
  return(SeuratObj)
}


#' Generate a marker file to run Garnett.
#'
#' @param markerDF A data frame with two columns (celltype, gene), containing the markers of each cell type.
#' @param SeuratObj The Seurat object that you going to annotate using Garnett.
#' @param file Path to the marker file to be generated. Default is "./markers.txt".
#'
#' @return None, marker file is written to \code{file} parameter location.
#' @export
#'
get_marker_file <- function(markerDF, SeuratObj, file = "./markers.txt") {
  markerDF <- markerDF[markerDF$gene %in% rownames(SeuratObj), ]
  sink(file)
  for (x in unique(markerDF$celltype)) {
    name <- gsub("\\(|\\)|:|>|,|#", ".", x)
    cat(">", name, "\n", sep = "")
    cat("expressed: ", paste(unique(markerDF$gene[markerDF$celltype == x]), collapse = ", "), "\n", sep = "")
    cat("\n")
  }
  sink()
}



#' Automatic annotation of Seurat object using scCATCH package.
#'
#' @param SeuratObj Seurat object to annotate
#' @param marker_custom Please refer to \url{https://raw.githack.com/ZJUFanLab/scCATCH/master/vignettes/tutorial.html} to
#' build a marker data frame. Or use \code{get_marker_scCATCH()} to generate the data frame.
#'
#' @return data frame of prediction result
#' @export
#'
AutoAnnotate_scCATCH <- function(SeuratObj, marker_custom) {

  obj <- scCATCH::createscCATCH(data = Seurat::GetAssayData(SeuratObj, slot = "data") ,
                                cluster = as.character(Seurat::Idents(SeuratObj)))
  obj <- scCATCH::findmarkergene(object = obj, if_use_custom_marker = TRUE, marker = marker_custom, use_method = "2")
  obj <- scCATCH::findcelltype(object = obj)
  result <- obj@celltype
  return(result)
}


#' Generate a data frame of marker genes with the required format to run scCATCH.
#'
#' @param markerDF A data frame with two columns (celltype, gene), containing the markers of each cell type.
#' @param species species.
#' @param tissue tissue.
#'
#' @return A data frame containing the markers of each cell type.
#' @export
#'
get_marker_scCATCH <- function(markerDF, species, tissue) {
  markerDF$species <- species
  markerDF$tissue <- tissue
  markerDF <- markerDF[, c('species', 'tissue', 'celltype', 'gene')]
  markerDF$pmid <- 11111111
  markerDF$subtype1 <- NA
  markerDF$subtype2 <- NA
  markerDF$subtype3 <- NA
  return(markerDF)
}









