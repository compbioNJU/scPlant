

#' visualize the Seurat object in a built-in shiny app
#'
#' @param SeuratObj Seurat object
#'
#' @export
#'
#' @examples
#' \dontrun{
#' if(interactive()) {
#'     load_shinyApp(SeuratObj)
#' }
#' }
#'
load_shinyApp <- function(SeuratObj, rasMat, rssMat, tf_target) {
  
  if (!inherits(x = SeuratObj, what = 'Seurat')) {
    stop("Seurat object is required, please provide a Seurat object.")
  }
  
  pkgs <- c('shiny', 'shinydashboard', 'data.table', 'DT', 'ggplot2', 'shinycssloaders', 'shinydashboardPlus', 'shinyWidgets',
            'dplyr', 'plyr', 'magrittr', 'Seurat', 'RColorBrewer', 'pheatmap', 'igraph', 'ggraph', 'circlize', 'ggpubr', 'cowplot',
            'factoextra', 'scatterpie', 'GO.db', 'grid', 'GOSemSim', 'ggplotify', 'GOfuncR', 'topicmodels', 'tidytext',
            'ggwordcloud', 'wordcloud', 'wordcloud2', 'pals', 'networkD3', 'echarts4r', 'ggnewscale')
  check_pkgs(pkgs)
  
  appDir <- system.file('shiny', package = "scPlant")
  shiny::runApp(appDir)
}




# Check packages for shiny app
#
# @param pkgs packages that need to check for shiny app
#' @importFrom utils install.packages
#
check_pkgs <- function(pkgs) {
  unavail <- c()
  for (pkg in pkgs) {
    if(!requireNamespace(pkg, quietly = TRUE)) {
      unavail <- c(unavail, pkg)
    }
  }
  if (length(unavail) > 0) {
    if(!interactive()) {
      stop("These packages are required but not installed, you need to manually install: \n",
           paste(strwrap(paste(unavail, collapse = ", ")), collapse = "\n"))
    }
    message("These packages are required but not installed: \n", paste(strwrap(paste(unavail, collapse = ", ")), collapse = "\n"))
    answer = readline("Do you want to install them? [y|n] ")
    
    if(tolower(answer) %in% c("y", "yes")) {
      for (pkg in unavail) {
        tryCatch(install.packages(pkg),
                 error = function(e) {
                   message("Installation from CRAN failed: \n", e$message, "\nTry to install ", pkg, " from Bioconductor")
                   if(!requireNamespace("BiocManager", quietly = TRUE)) {
                     install.packages("BiocManager")
                   }
                   BiocManager::install(pkg)
                 })
      }
    } else {
      stop("You need to manually install these packages.")
    }
  }
}


