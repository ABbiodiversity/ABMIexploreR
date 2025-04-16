#' Coefficient extraction
#'
#' @description This function extracts the stored coefficients and provides a table or summary figure.
#' 
#' @param species Unique Species ID. If defined as "All", the complete coefficient table will be provided to the user.
#' @param model Vector defined as "vegetation" or "soil" that provides the relevant coefficients.
#'
#' @import Matrix
#'
#' @export
#'
#'

coefficient_extraction <- function(species, model = NULL) {
    
    # Check for valid model
    if(!(model %in% c("vegetation", "soil"))) {
        
        stop("Need to a valid model type (vegetation, soil).")
    }
    
    # Generate the species lookup table
    species.lookup <- abmi_species()
    
    if(model == "vegetation") {
        
        species.lookup <- species.lookup[species.lookup$ModelNorth, ]
        
    }
    
    if(model == "soil") {
        
        species.lookup <- species.lookup[species.lookup$ModelSouth, ]
        
    }
    
    # Generate the lookup table
    coef.table <- NULL
    for(species.id in species.lookup$SpeciesID) {
        
        # Isolate the relevant bootstrap iteration
        boot <- species.lookup[species.lookup$SpeciesID %in% species.id, "Bootstrap"]
        taxon <- species.lookup[species.lookup$SpeciesID %in% species.id, "Taxon"]
        
        # Extract the coefficient table
        coef.template <- species.coefs[[taxon]][[model]][species.id, , boot]
        coef.template <- coef.template[!(names(coef.template) %in% c("Climate", "pAspen"))]

        # Backtransform coefficients based on model type
        inv_link <- inv_link_function[[taxon]]
        coef.template <- inv_link(coef.template)
        
        # Format
        coef.template <- matrix(coef.template, nrow = 1, dimnames = list(species.id, names(coef.template)))
        
        coef.table <- rbind(coef.table, coef.template)
        
    }
    
    # If all species, generate the entire lookup table
    if(species == "All") {
        
        return(coef.table)
        
    } 
    
    # If only a single species, generate the single set of coefficients
    if(species %in% species.lookup$SpeciesID) {
        
        coef.table <- coef.table[species, ]
        return(coef.table)
        
    } else {
        
        stop("Species model not available. Check SpeciesID or model type (vegetation or soil).")
        
    }
    
}

