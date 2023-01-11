source(
  file = "global.R",
  local = TRUE,
  encoding = "UTF-8"
)

server <- function(input, output, session) {
  
  
  source(
    file = "function/function.server.R",
    local = TRUE,
    encoding = "UTF-8"
  )
  
  source(
    file = "DEG/DEG.server.R",
    local = TRUE,
    encoding = "UTF-8"
  )

  source(
    file = "network/network.server.R",
    local = TRUE,
    encoding = "UTF-8"
  )
  
  


}





