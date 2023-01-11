fluidPage(
  fluidRow(
    div(class="customdiv", 
        h2(align="left", "GSEA result"),
      fluidRow(
        column(
          width = 11,
          radioGroupButtons(
            justified = F,
            inputId = "switchTerm",
            label = NULL,
            # choices = c("GO", "KEGG", "Reactome", "MSigDb", "WikiPathways", "DO", "NCG", "DGN", "customized"),
            choices = c("GO"),
            # selected = "GO",
            checkIcon = list(
              yes = icon("ok", lib = "glyphicon")
            ),
            status = "primary",
            individual = T
          )
        ),
        column(
          width = 1,
          uiOutput("grouping")
        )
      ),
      DT::dataTableOutput("GSEAtable") %>% withSpinner(),
      # uiOutput("select_plotui") %>% withSpinner()
      fluidRow(
        id = "tohide2",
        column(width = 12,
               uiOutput("select_plotui2") %>% withSpinner()
        )
      ),
      fluidRow(
        id = "tohide1",
        br(),
        column(width = 4,
               uiOutput("pathwayList")
        ),
        column(width = 8,
               actionBttn(
                 inputId = "button2",
                 label = "Clear Selection",
                 style = "simple",
                 color = "primary",
                 icon=icon("broom"),
                 size = "sm"
               ),
               actionBttn(
                 inputId = "button4",
                 label = "Visualize single pathway",
                 icon=icon("arrow-circle-right"),
                 style = "unite",
                 color = "warning",
                 size = "sm"
               ),
               actionBttn(
                 inputId = "button3",
                 label = "Visualize multiple pathways",
                 style = "jelly",
                 icon=icon("arrow-circle-right"),
                 color = "danger",
                 size = "sm"
               )
        )
      )
    )
  ),
  br(),
  fluidRow(
    id = "pa",
    column(
      width = 6,
      div(class="customdiv", 
          gseaheatmapui(),
          plotOutput('gseaheatmap', height = "950px") %>% withSpinner()
      )
    ),
    column(
      width = 6,
      div(class="customdiv", 
          pathwayScatterplotui(),
          plotOutput('pathwayScatterplot', height = "450px") %>% withSpinner()
      ),
      div(class="customdiv", 
          edgeBundlingui(),
          plotOutput('edgeBundling', height = "450px") %>% withSpinner()
      )
    )
    # br(),
    # fluidRow(
    #   column(
    #     width = 6,
    #     div(class="customdiv", 
    #         h2(align="left", "Cosine_networkByGSEA"),
    #         cosineByGSEAui(),
    #         plotOutput('cosineByGSEA', height = "600px") %>% withSpinner()
    #     )
    #   ),
    #   column(
    #     width = 6,
    #     div(class="customdiv", 
    #         h2(align="left", "Clustercorplot"),
    #         clustercorplotui(),
    #         plotOutput('clustercorplot', height = "600px") %>% withSpinner()
    #     )
    #   )
    # )
  ),
  br(),
  fluidRow(
    column(width = 12,
           div(class="customdiv", 
               hierarchyplot_treeui(),
               plotOutput('hierarchyplot_tree', height = "900px") %>% withSpinner()
           )
    )
  ),
  br(),
  fluidRow(
    div(class="customdiv", 
        h2(align="left", tagList(icon("project-diagram"), "Topic modeling network")),
        hr(),
        # p(tagList(icon("lightbulb"), "Click network to view umap and terms of specific topic.")),
        fluidRow(
          column(width = 8,
                 numericInput(inputId = "tm_topn1", "show top terms:", 8, min = 1, max = 50, width = "300px"),
                 echarts4rOutput("tm_echt", height = "1000px") %>% withSpinner()
          ),
          column(width = 4,
                 div(style="display: flex;flex-direction: row;",
                   selectInput(
                     inputId = "tm_topics",
                     label = "Show Topic:",
                     width = '200px',
                     choices = 1:34
                   ),
                   div(style="display: flex;align-items: center;padding: 5px;",
                     actionBttn(
                       inputId = "submit",
                       label = "submit",
                       style = "simple",
                       color = "primary",
                       size = "sm"
                     )
                   )
                 ),
                 plotOutput('tm_topicProb') %>% withSpinner(),
                 br(),
                 plotOutput('tm_Topterms') %>% withSpinner()
          )
        )
    )
  ),
  br(),
  fluidRow(
    div(class="customdiv", 
        fluidRow(
          column(width = 5,
                 sankeyNetworkOutput("tm_sankey", height = "1200px") %>% withSpinner()
          ),
          column(width = 7,
                 plotOutput('tm_heatmap', height = "600px") %>% withSpinner(),
                 plotOutput('tm_cosine_network_cluster', height = "550px") %>% withSpinner(),
                 fluidRow(
                   column(width = 6,
                          sliderInput("cos_sim_thresh", label = "cosine similarity threshold", value = 0.2,
                                      min = 0.01, max = 1, step = 0.01)
                   ),
                   column(width = 6,
                          sliderInput("radiuscosnet", label = "pie radius", value = 0.15,
                                      min = 0.01, max = 1, step = 0.01)
                   )
                 )
          )
        )
    )
  ),
  br(),
  fluidRow(
    div(class="customdiv", 
        fluidRow(
          column(width = 2,
                 div(style="padding: 10px;",
                   sliderInput("tm_cosnet_topn", label = "show top terms of every topic/cluster", value = 10,
                               min = 1, max = 50, step = 1),
                   selectInput(
                     inputId = "tm_layout2",
                     label = "layout", 
                     selected = "fr",
                     choices = c("fr", "kk", "star", "circle", "dh", "graphopt", "graphopt", "sphere", "drl", "nicely", "lgl")
                   ),
                   sliderInput("cos_sim_thresh2", label = "cosine similarity threshold", value = 0.8,
                               min = 0.01, max = 1, step = 0.01),
                   sliderInput("radiuscosnet2", label = "pie radius", value = 0.2,
                               min = 0.01, max = 1, step = 0.01),
                   sliderInput("tm_text_size", label = "text size", value = 3,
                               min = 0.5, max = 10, step = 0.5)
                 )
          ),
          column(width = 10,
                 plotOutput('tm_cosine_network_term', height = "1000px") %>% withSpinner()
          )
        )
    )
  )
  
)






