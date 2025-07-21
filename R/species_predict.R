#' Species predictions
#'
#' @description Wrapper function for calculating the species level predictions
#'
#' @param species Unique Species ID. 
#' @param veg Matrix of standardized vegetation information (row sum equals 1). If NULL and model is available, no prediction is generated.
#' @param soil Matrix of standardized soil information (row sum equals 1). If NULL and model is available, no prediction is generated.
#' @param climate Matrix of the bioclimatic variables that is matched to the vegetation and/or soil matrices.
#' @param modified Logical; Default behaviour (FALSE) means the original species coefficients are used. TRUE indicators the function should use the coefficients created through the coefficient_adjustment function.
#' @param boot Defines bootstrap iteration to generate the predictions from. If 0, uses the median bootstrap run unique for each species.
#'
#' @import Matrix
#'
#' @export
#'

species_predict <- function(species, veg = NULL, soil = NULL, 
                            climate = NULL, modified = FALSE, boot = 0) {
    
    ##########################
    # Check for valid inputs #
    ##########################
    
    if(is.null(veg) & is.null(soil)) {
        
        stop("Need to define the landcover information.")
        
    }
    
    if(is.null(climate)) {
        
        stop("Need to define the climate information.")
        
    }
    
    # Check that there are no mismatches between the landcover and coefficients
    if(!is.null(veg)) {
        
        landcover.lookup <- landcover.names$Vegetation
        
        if(!all(colnames(veg) %in% landcover.lookup$Variable[landcover.lookup$Type == "Vegetation"])) {
            
            stop(paste("Unknown vegetation or footprint classes detected: %s",
                       paste(setdiff(colnames(veg), landcover.lookup$Variable[landcover.lookup$Type == "Vegetation"]),
                             collapse=", ")))
            
        }
        
    }
    
    if(!is.null(soil)) {
        
        landcover.lookup <- landcover.names$Soil
        
        if(!all(colnames(soil) %in% landcover.lookup$Variable[landcover.lookup$Type == "Soil"])) {
            
            stop(paste("Unknown soil or footprint classes detected: %s",
                       paste(setdiff(colnames(soil), landcover.lookup$Variable[landcover.lookup$Type == "Soil"]),
                             collapse=", ")))
            
        }
        
    }
    
    ###############################
    # Define available bootstraps #
    ###############################
    
    taxon <- .get_taxon(species)
    
    if (boot > 100 | boot < 0) {
        
        .msg("Bootstrap ID must be between 0-100", 4)
        
    }
    
    if (taxon %in% c("Mammals") && boot > 1) {
        
        .msg(sprintf("Bootstrap coefs (boot > 1) not available for species %s (%s)", species, taxon), 4)
        
    } 
    
    #####################################
    # Make the appropriate climate data #
    #####################################
    
    taxon.climate <- .create_climate(climate.grid = climate, taxon = taxon) 
    
    ###############################
    # Align veg with climate data #
    ###############################
    
    veg <- veg[rownames(taxon.climate$vegetation), ]
    soil <- soil[rownames(taxon.climate$soil), ]
    
    ###########################
    # Define prediction model #
    ###########################
    
    .msg(sprintf("Making predictions for species %s (%s) boot=%s", species, taxon, boot))
    
    if (taxon == "Mammals") {
        
        species.pred <- mammal_predicts(species = species, taxon = taxon, veg = veg, 
                                        soil = soil, climate = taxon.climate, 
                                        boot = boot, modified = modified)
        
    }
    
    if (taxon %in% c("Amphibians", "Birds", "Bryophytes", "Lichens", "Mites", "VascularPlants")) {
        
        species.pred <- .terrestrial_predicts(species = species, taxon = taxon, veg = veg, 
                                              soil = soil, climate = taxon.climate, 
                                              boot = boot, modified = modified)
        
    }
    
    # Return the prediction
    return(species.pred)
    
}

