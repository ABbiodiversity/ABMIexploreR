#' Internal lookup tables
#'
#' @section Provide users with access to internal lookup tables.
#'
#' @description This function provides the lookup table of the vegetation, soil, and footprint definitions
#'
#' @param landcover Define 'Vegetation' or 'Soil' landcover classes.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' landcover.lookup <- abmi_landcover(landcover = "Vegetation")
#' }
#'
#' @return A dataframe with all of the coefficients used in the joint climate and landcover model.
#' 

abmi_landcover <- function(landcover) {
    
    return(landcover.names[[landcover]])
    
}

#' Provide users with bioclimatic lookup table
#'
#' @description This function provides the lookup table of the bioclimatic covariates used in the climate model.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' climate.lookup <- abmi_climate()
#' }
#'
#' @return A dataframe with all of the coefficients used in the climate model. Each taxonomic group uses a unique set of these variables.

abmi_climate <- function() {
    
    return(bioclimatic.lookup)
    
}

#' Provide users with species lookup table
#'
#' @description This function provides the lookup table of all species models.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' species.lookup <- abmi_climate()
#' }
#'
#' @return A dataframe with all of the species information include model region and validation statistics.

abmi_species <- function() {
    
    species.lookup <- NULL
    for(taxon in names(species.coefs)) {
        
        species.lookup <- rbind(species.lookup, species.coefs[[taxon]]$Species)
        
    }
    
    return(species.lookup)
    
}
