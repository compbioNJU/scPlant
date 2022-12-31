

get_orgDb <- function(species = c('Ath', 'Osa', 'Zma')) {
  
  species <- match.arg(species)
  if (species == 'Ath') {
    suppressPackageStartupMessages(require('org.At.tair.db', character.only = TRUE))
    orgdb <- org.At.tair.db
  } else if (species == 'Osa') {
    orgdb <- AnnotationDbi::loadDb(file = system.file("extdata", "Oryza_sativa_orgdb.sqlite", package = "scPlant"))
  } else {
    orgdb <- AnnotationDbi::loadDb(file = system.file("extdata", "Zea_mays_orgdb.sqlite", package = "scPlant"))
  }
  return(orgdb)
}


