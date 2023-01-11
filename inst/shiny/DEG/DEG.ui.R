fluidPage(
  fluidRow(
    box(width = 12,
        title = "Differentially Expressed Genes",
        br(),
        fluidRow(
          column(
            width = 12,
            uiOutput("geneTable") %>% withSpinner(),
          )
        ),
        br(),
        fluidRow(
          column(4,
                 uiOutput("geneList")
          ),
          column(8,
                 actionBttn(
                   inputId = "btn1",
                   label = "Clear Selection",
                   style = "simple",
                   color = "primary",
                   icon=icon("broom"),
                   size = "sm"
                 ),
                 actionBttn(
                   inputId = "btn2",
                   label = "Show single gene",
                   icon=icon("arrow-circle-right"),
                   style = "unite",
                   color = "warning",
                   size = "sm"
                 ),
                 actionBttn(
                   inputId = "btn3",
                   label = "Show multiple genes",
                   icon=icon("arrow-circle-right"),
                   style = "jelly",
                   color = "danger",
                   size = "sm"
                 )
          )
        )
      )
    
  ),
  br(),
  fluidRow(
    box(width = 12,
      tabsetPanel(
        id = "tabset",
        tabPanel(title =strong("single gene"), value = "single",
                 fluidRow(
                   column(width = 3,
                          box(width = 12,
                              plotOutput("singleGene_featureplot", width = "450px", height="400px") %>% withSpinner()
                          )
                   ),
                   column(width = 9,
                          box(width = 12,
                              plotOutput("singleGene_violinplot", height="400px") %>% withSpinner()
                          )
                   )
                 )
        ),
        tabPanel(title = strong("multiple genes"), value = "multiple",
                 fluidRow(
                   box(width = 12,
                       h3(align="center", "Bubble Plot"),
                       plotOutput("bubble_plot", height="500px") %>% withSpinner()
                   ),
                   box(width = 12,
                       h3(align="center", "Violin Plot"),
                       plotOutput("vln_plot", height="500px") %>% withSpinner()
                   )#,
                   # box(width = 12,
                   #     h3(align="center", "Feature Plot"),
                   #     plotOutput("feature_plot", height="500px") %>% withSpinner()
                   # )
                 )
        )
      )
    )

  )
)



