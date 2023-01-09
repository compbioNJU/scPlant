

#' Decompose bulk RNA-seq samples using a single cell reference
#'
#' CIBERSORT is an algorithm (Newman et al. Nature Methods 12:453-457) for estimating the cell type composition of a bulk sample,
#' given a gene expression profile of the sample and a known gene expression profile for each cell type potentially contributing to the sample.
#' Here, we wrap runCIBERSORT function from RNAMagnet R package to decompose bulk RNA-seq samples.
#' To install RNAMagnet R package, see https://github.com/veltenlab/rnamagnet. Note that MAGIC is an important requirement of RNAMagnet package,
#' and needs to be installed as a python package, see https://github.com/KrishnaswamyLab/MAGIC/tree/master/Rmagic#installation.
#'
#' @param SeuratObj A Seurat object.
#' @param exprs A data frame or matrix of raw read counts of bulk RNA-seq samples. Column names correspond to sample names, row names to genes.
#' @param meta A data frame of 2 columns with Sample and Class, correlate sample names with sample class.
#' @param mc.cores Number of cores used
#' @param slot Which slot to use to calculate mean expression per cell type
#' @param assay Which assays to use to calculate mean expression per cell type
#'
#' @return A data frame in long format, indicating the fraction of cell types in each bulk RNA-seq sample.
#' @export
#'
decomposeBulk <- function(SeuratObj, exprs, meta, mc.cores = 1, slot = "counts", assay = 'SCT') {
  # 1. identify marker genes from a seurat object. Here we chose top 20 markers.
  if (!"Allmarkers" %in% names(slot(object = SeuratObj, name = 'misc'))) {
    message("Identifying marker genes from the single cell reference......")
    SeuratObj.markers <- Seurat::FindAllMarkers(SeuratObj, only.pos = TRUE, min.pct = 0.0001, logfc.threshold = 0.0001, return.thresh=0.9)
    slot(object = SeuratObj, name = 'misc')[["Allmarkers"]] <- SeuratObj.markers
  }
  Allmarkers <- slot(object = SeuratObj, name = 'misc')[["Allmarkers"]]
  usegenes <- Allmarkers %>% dplyr::group_by(cluster) %>% dplyr::slice_max(order_by = avg_log2FC, n = 20, with_ties = F) %>%
    dplyr::pull(gene) %>% unique
  # 2. compute mean expression per cell type
  mean_by_cluster <- Seurat::AverageExpression(SeuratObj, features = usegenes, slot = slot, assay = assay)[[assay]]
  mean_by_cluster <- as.matrix(mean_by_cluster)
  # 3. Create a vector that maps samples to biological class
  meta <- meta[, c("Sample", "Class")] %>% tibble::deframe()
  meta <- meta[colnames(exprs)]
  # 4. Run CIBERSORT
  CIBER <- runCIBERSORT(exprs=exprs, base=mean_by_cluster, markergenes=intersect(rownames(exprs), rownames(mean_by_cluster)),
                        design=meta, mc.cores=mc.cores)
  return(CIBER)
}

#' modified version of original code(https://rdrr.io/github/veltenlab/rnamagnet/src/R/CIBERSORT.R).
runCIBERSORT <- function(exprs, base, design, markergenes = intersect(rownames(base), rownames(exprs)),
                         transform=function(x) x,nu = c(0.25,0.5,0.75), optim.nu = F, mc.cores= 3, ...) {

  res <- list()
  for (i in 1:ncol(exprs)) {
    x <- exprs[,i]
    names(x) <- rownames(exprs)
    res[[i]] <- RNAMagnet::CIBERSORT(x, features=base, transform=transform, usegenes = intersect(markergenes, rownames(exprs)), nu=nu, optim.nu = optim.nu, mc.cores = mc.cores, ...)
  }
  #out <- apply(exprs,2, CIBERSORT, features=mean_by_cluster, kernel=kernel, cost =cost, method = method,alpha=alpha, gamma=gamma, transform=transform, usegenes = intersect(markergenes, rownames(exprs)), norm=norm, nu=nu)
  pvals <- data.frame(
    pvals = sapply(res, attr, "p"),
    samples = colnames(exprs)
  )

  svs <- lapply(res, attr, "SV")

  out <- do.call(cbind,res)
  colnames(out) <- colnames(exprs)
  rownames(out) <- colnames(base)
  out <- reshape2::melt(out)
  out$experiment <- design[out$Var2]
  colnames(out) <- c("CellType","SampleID","Fraction","SampleClass")
  out
}




#' Heatmap showing the fraction of cell types in each bulk RNA-seq sample.
#'
#' @param CIBER output of \code{decomposeBulk()}, a data frame indicating the fraction of cell types in each bulk RNA-seq sample.
#' @param meta A data frame of 2 columns with Sample and Class, correlate sample names with sample class.
#'
#' @importFrom ComplexHeatmap Heatmap rowAnnotation
#' @importFrom grid gpar
#'
#' @export
#'
fractionHmp <- function(CIBER, meta) {
  meta <- meta[, c("Sample", "Class")] %>% tibble::deframe()
  mmm <- CIBER %>% dplyr::mutate(CellType=as.character(CellType)) %>% reshape2::acast(SampleID~CellType, value.var = "Fraction")
  col_fun = circlize::colorRamp2(seq(min(mmm), max(mmm), length.out = 10), c("white", pals::brewer.orrd(9)))
  ht <- Heatmap(mmm[names(meta), ], name = "Fraction", col = col_fun, column_names_rot = 45, row_names_gp = gpar(fontsize = 7),
                cluster_rows = F, cluster_columns = T,
                column_names_gp = gpar(fontsize = 7),  row_split=meta, row_title_rot = 0, row_title_gp = gpar(fontsize = 8),
                column_names_centered = F, border = TRUE, rect_gp = gpar(col = NA), show_row_names = F,
                left_annotation = rowAnnotation(sample = meta, width = unit(1, "mm")))
  return(ht)
}


#' Heatmap showing mean expression of top 20 markers of each cell type in each sample.
#'
#' @param SeuratObj A Seurat object.
#' @param exprs A data frame or matrix of raw read counts of bulk RNA-seq samples. Column names correspond to sample names, row names to genes.
#' @param meta A data frame of 2 columns with Sample and Class, correlate sample names with sample class.
#'
#' @importFrom ComplexHeatmap Heatmap rowAnnotation draw
#' @importFrom grid gpar
#'
#' @export
#'
meanExpHmp <- function(SeuratObj, exprs, meta) {
  meta <- meta[, c("Sample", "Class")] %>% tibble::deframe()
  exprs <- exprs[, names(meta)]
  exprs <- exprs[rowSums(exprs)>0, ]
  exprs <- t(apply(exprs, 1, function(x){(x-min(x))/(max(x)-min(x))}))

  Allmarkers <- slot(object = SeuratObj, name = 'misc')[["Allmarkers"]]
  topMarker <- Allmarkers %>% dplyr::group_by(cluster) %>% dplyr::slice_max(order_by = avg_log2FC, n = 20, with_ties = F) %>% dplyr::ungroup()
  meanExp <- sapply(split(topMarker$gene, topMarker$cluster), function(gg){
    gg <- gg[gg %in% rownames(exprs)]
    colMeans(exprs[unique(gg), ])
  })
  meanExp <- t(scale(t(meanExp), center = T, scale=T))

  col_fun = circlize::colorRamp2(seq(min(meanExp), max(meanExp), length.out = 10), c("white", pals::brewer.blues(9)))
  ht <- Heatmap(meanExp[names(meta), ], name = "mean expression", col = col_fun, column_names_rot = 45, row_names_gp = gpar(fontsize = 7),
                cluster_rows = F, cluster_columns = T,
                column_names_gp = gpar(fontsize = 7),
                row_split=meta, row_title_rot = 0, row_title_gp = gpar(fontsize = 8),
                column_names_centered = F, border = TRUE, rect_gp = gpar(col = NA), show_row_names = F,
                left_annotation = rowAnnotation(sample = meta, width = unit(1, "mm")))
  return(ht)
}





