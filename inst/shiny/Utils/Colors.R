
palettes <- function() {
  return(list("Greens", "Blues", "BuGn", "BuPu", "GnBu", "Greys", "Oranges", "OrRd", "PuBu",
              "PuBuGn", "PuRd", "Purples", "RdPu", "Reds", "YlGn", "YlGnBu", "YlOrBr", "YlOrRd"))
}

color_choicesOpt <- function(palettes, id) {
  switch(id,
    `1` = list(content = sapply(palettes,
                              function(x){
                                paste0("<div style='background: ", brewer.pal(7,x)[7], "; color: white; padding-left: 5px; padding-right: 5px;'>", x, "</div>")
                              })),
    `2` = list(style = sapply(palettes,
                            function(x){
                              paste0("background: ", brewer.pal(7,x)[7], "; color: white;")
                            })),
    `3` = list(style = sapply(palettes,
                            function(x){
                              paste0("color: black; font-weight: bold; background-image: linear-gradient(to right, ",
                                     paste(c('white', brewer.pal(n=7,name=x)), collapse = ", "), ");")
                              # paste0("color: black; font-weight: bold; background-image: linear-gradient(to right, ",
                              #        paste(colorRampPalette(c('white', brewer.pal(n=7,name=x)))(100), collapse = ", "), ");")
                            })),
    `4` = list(style = sapply(palettes,
                            function(x){
                              paste0("color: ", brewer.pal(7,x)[7], "; font-weight: bold;")
                            }))
  )
}




