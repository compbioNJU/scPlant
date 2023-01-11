fluidPage(
  fluidRow(
    box(
      title = "Top regulons of each celltype",
      solidHeader = T,
      width = 12,
        fluidRow(
          column(
            width = 4,
            DT::dataTableOutput("RSStable") %>% withSpinner(),
            br(),
            fluidRow(
              column(width = 12,
                     strong("Select rows to:"),
                     actionBttn(
                       inputId = "nw_btn1",
                       label = "show these regulons",
                       icon=icon("arrow-circle-right"),
                       style = "unite",
                       color = "warning",
                       size = "sm"
                     ),
                     actionBttn(
                       inputId = "nw_btn2",
                       label = "Clear Selection",
                       style = "simple",
                       color = "primary",
                       icon=icon("broom"),
                       size = "sm"
                     )
                     # actionBttn(
                     #   inputId = "nw_btn3",
                     #   label = "View this regulon",
                     #   style = "jelly",
                     #   icon=icon("arrow-circle-right"),
                     #   color = "danger",
                     #   size = "sm"
                     # )
              )
            )
          ),
          column(
            width = 8,
            numericInput(inputId = "nw_topn1", "show top regulons of each celltype:", 10, min = 1, max = 50, width = "500px"),
            plotOutput('RSS_topnPlot', height = "1000px") %>% withSpinner()
          )
        )
    )
  ),
  # fluidRow(
  #   id = "nw_tohide2",
  #   box(
  #       hr(),
  #       fluidRow(
  #         column(
  #           width = 4,
  #           tags$img(src = "umap.png", width = "400px", height = "400px")
  #         ),
  #         column(
  #           width = 8,
  #           plotOutput('ras_exp_plot', height = "400px") %>% withSpinner()
  #         )
  #       )
  #   )
  # ),
  hr(),
  fluidRow(
    id = "nw_tohide1",
    box(
      title = "Top targets of each regulon",
      width = 12,
        fluidRow(
          column(
            width = 4,
            DT::dataTableOutput("targetTable") %>% withSpinner()
          ),
          column(
            width = 8,
            numericInput(inputId = "nw_topn2", "show top targets of each regulon:", 5, min = 1, max = 20, width = "500px"),
            plotOutput('target_topnPlot', height = "700px") %>% withSpinner()
          )
        ),
        hr(),
        fluidRow(
          id = "nw_tohide3",
          column(width = 12,
                 plotOutput('ras_exp_hmpPlot', height = "600px") %>% withSpinner()
                 )
        )
    )
  )
)












