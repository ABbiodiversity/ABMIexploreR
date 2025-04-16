#' Determine taxonomic group
#'
#' @description Function for determining a species taxonomic group
#'
#' @param species Valid species name available from the abmi_species function.
#' @param taxon Taxonomic group associated with the focal species of interest. 
#'
#' @export
#'

.get_taxon <- function(species) {

    # Call the species lookup table function
    species.lookup <- abmi_species()
    
    # Isolate the taxonomic group
    taxon <- species.lookup[species.lookup$SpeciesID %in% species, "Taxon"]
    
    if(length(taxon) == 0) {
        
        stop(paste0("Coefs not available for species: ", species))
        
    } else {
        
        return(taxon)
        
    }
    

}
