
# clustercorplotui <- function() {
#   fluidRow(
#     column(width = 10),
#     column(width = 2,
#            tags$span(
#              dropdownButton(
#                tags$h4(icon("sliders-h"), "Options"),
#                tags$hr(),
#                radioGroupButtons(
#                  inputId = "similarity",
#                  label = "Similarity measurement",
#                  choices = c("Jaccard index", "Pearson's correlation"),
#                  status = "primary",
#                  checkIcon = list(
#                    yes = icon("check-square"),
#                    no = icon("square-o")
#                  ),
#                  selected = "Pearson's correlation",
#                  individual = F,
#                  direction = "vertical"
#                ),
#                sliderInput("link_threshold_clustercorplot", label = "only show links whose correlation/Jaccard-index bigger than", value = 0.5,
#                            min = 0, max = 1, step = 0.05),
#                sliderInput("vertex.size.cex", label = "vertex.size.cex", value = 1,
#                            min = 0, max = 5, step = 0.1),
#                sliderInput("vertex.label.cex", label = "vertex.label.cex", value = 1,
#                            min = 0, max = 5, step = 0.1),
#                sliderInput("edge.max.width", label = "edge.max.width", value = 8,
#                            min = 1, max = 30, step = 1),
#                sliderInput("vertex.label.dist", label = "vertex.label.dist", value = 2,
#                            min = 0, max = 10, step = 0.5),
#                downloadButton("clustercorplotdl", label="Download PDF"),
#                circle = F,
#                right = T,
#                inline = T,
#                status = "primary",
#                icon = icon("gear"),
#                size = "sm",
#                width = "300px",
#                tooltip = tooltipOptions(title = "Options", placement = "top")
#              ),
#              style = "float:right;"
#            )
#     )
#   )
# }


gseaheatmapui <- function() {
  fluidRow(
    column(width = 7),
    column(width = 4,
           pickerInput(
             inputId = "hmpcol",
             label = "color",
             choices = palettes(),
             choicesOpt = color_choicesOpt(palettes(), id=3),
             inline = T
           )
    ),
    column(width = 1,
           tags$span(
             dropdownButton(
               tags$h4(icon("sliders-h"), "Options"),
               tags$hr(),
               sliderInput("topPath", label = "top pathways of every cluster to show", value = 6,
                           min = 1, max = 30, step = 1),
               sliderInput("fontsize_row", label = "row fontsize", value = 10,
                           min = 1, max = 30, step = 1),
               selectInput("toshow", label = "heatmap showing",
                           choices = list("-log10(p.adjust)" = "-logFDR", "normalized enrichment score(NES)" = "NES"),
                           selected = "-logFDR"),
               radioGroupButtons(
                 inputId = "scale",
                 label = "scale",
                 choices = c("row", "none", "column"),
                 status = "primary",
                 individual = T
               ),
               downloadButton("gseaHeatmapdl", label="Download PDF"),
               # sliderInput("plotheight", label = "height of plot", value = 800,
               #             min = 100, max = 2000, step = 50),
               circle = F,
               right = T,
               inline = T,
               status = "primary",
               icon = icon("gear"),
               size = "sm",
               width = "300px",
               tooltip = tooltipOptions(title = "Options", placement = "top")
             ),
             style = "float:right;"
           )
    )
  )
}





hierarchyplot_treeui <- function() {
  fluidRow(
    column(width = 10),
    column(width = 2,
           tags$span(
             dropdownButton(
               tags$h4(icon("sliders-h"), "Options"),
               tags$hr(),
               sliderInput("topath1", label = "top pathways of every cluster to show", value = 3,
                           min = 1, max = 50, step = 1),
               numericInput(inputId = "cluster_cutree_k", label = "cluster_cutree_k", 14, min = 1, max = 50),
               numericInput(inputId = "pathway_cutree_k", label = "pathway_cutree_k", 14, min = 1, max = 50),
               sliderInput("edge.max.width2", label = "edge.max.width", value = 1.5,
                           min = 0, max = 5, step = 0.5),
               sliderInput("vertex.size.cex2", label = "vertex.size.cex", value = 0.5,
                           min = 0, max = 5, step = 0.1),
               sliderInput("vertex.label.cex2", label = "vertex.label.cex", value = 0.8,
                           min = 0, max = 5, step = 0.1),
               sliderInput("alpha.edge", label = "transparency of edge", value = 0.6,
                           min = 0, max = 1, step = 0.05),
               downloadButton("hierarchyplot_treedl", label="Download PDF"),
               circle = F,
               right = T,
               inline = T,
               status = "primary",
               icon = icon("gear"),
               size = "sm",
               width = "300px",
               tooltip = tooltipOptions(title = "Options", placement = "top")
             ),
             style = "float:right;"
           )
    )
  )
}




pathwayScatterplotui <- function() {
  fluidRow(
    column(width = 7),
    column(width = 4,
           pickerInput(
             inputId = "sctcol",
             label = "color",
             choices = palettes(),
             selected = "OrRd",
             choicesOpt = color_choicesOpt(palettes(), id=3),
             inline = T
           )
    ),
    column(width = 1,
           tags$span(
             dropdownButton(
               tags$h4(icon("sliders-h"), "Options"),
               tags$hr(),
               materialSwitch(
                 inputId = "labelornot",
                 label = tags$b("Label the clusters"),
                 status = "primary",
                 value = TRUE
               ),
               sliderInput("pointsize", label = "point size", value = 1,
                           min = 0, max = 10, step = 0.1),
               sliderInput("label.size", label = "label size", value = 4,
                           min = 1, max = 10, step = 0.2),
               radioGroupButtons(
                 inputId = "reduction1",
                 label = "reduction",
                 choices = c("umap", "tsne", "pca"),
                 status = "primary",
                 individual = T
               ),
               downloadButton("pathwayScatterplotdl", label="Download PDF"),
               circle = F,
               right = T,
               inline = T,
               status = "primary",
               icon = icon("gear"),
               size = "sm",
               width = "300px",
               tooltip = tooltipOptions(title = "Options", placement = "top")
             ),
             style = "float:right;"
           )
    )
  )
}


# cosineByGSEAui <- function() {
#   fluidRow(
#     column(width = 10),
#     column(width = 2,
#            tags$span(
#              dropdownButton(
#                tags$h4(icon("sliders-h"), "Options"),
#                tags$hr(),
#                sliderInput("cosineByGSEA_thr", label = "only show links with cosine similarity bigger than", value = 0.7,
#                            min = 0, max = 1, step = 0.05),
#                sliderInput("cosineByGSEA_text", label = "text size", value = 0.8,
#                            min = 0, max = 3, step = 0.1),
#                circle = F,
#                right = T,
#                inline = T,
#                status = "primary",
#                icon = icon("gear"),
#                size = "sm",
#                width = "300px",
#                tooltip = tooltipOptions(title = "Options", placement = "top")
#              ),
#              style = "float:right;"
#            )
#     )
#   )
# }


edgeBundlingui <- function() {
  fluidRow(
    column(width = 10),
    column(width = 2,
           tags$span(
             dropdownButton(
               tags$h4(icon("sliders-h"), "Options"),
               tags$hr(),
               radioGroupButtons(
                 inputId = "edgeBundling_method",
                 label = "links showing:",
                 choices = c('cosine similarity', 'pearson', 'spearman'),
                 status = "primary",
                 checkIcon = list(
                   yes = icon("check-square"),
                   no = icon("square-o")
                 ),
                 selected = "cosine similarity",
                 individual = F,
                 direction = "vertical"
               ),
               sliderInput("edgeBundling_thr", label = "only show links with similarity/correlation bigger than", value = 0.7,
                           min = 0, max = 1, step = 0.05),
               circle = F,
               right = T,
               inline = T,
               status = "primary",
               icon = icon("gear"),
               size = "sm",
               width = "300px",
               tooltip = tooltipOptions(title = "Options", placement = "top")
             ),
             style = "float:right;"
           )
    )
  )
}







