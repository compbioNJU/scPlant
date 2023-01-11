source(
  file = "global.R",
  local = TRUE,
  encoding = "UTF-8"
)

ui <- navbarPage(title = "scPlant-App",
                 header=tagList(
                   shinyjs::useShinyjs(),
                   singleton(tags$head(tags$link(rel="stylesheet", type = "text/css", href = "custom.css")))
                   ),
                 theme = shinytheme("cerulean"),
                 
                 tabPanel(title = "Functional Annotation",
                          source(
                            file = "function/funciton.ui.R",
                            local = TRUE,
                            encoding = "UTF-8"
                          )$value
                          ),
                 tabPanel(title = "Differential Expression",
                          source(
                            file = "DEG/DEG.ui.R",
                            local = TRUE,
                            encoding = "UTF-8"
                          )$value
                 ),
                 tabPanel(title = "Regulatory Network",
                          source(
                            file = "network/network.ui.R",
                            local = TRUE,
                            encoding = "UTF-8"
                          )$value
                 )
                 
)
