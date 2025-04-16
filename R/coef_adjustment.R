#' Coefficient Adjustment
#'
#' @description This function will adjust the stored coefficients to user defined values. 
#' In general, we advice not adjusting these coefficients. However, there are specific situations 
#' where adjusting stored coefficients to 0 (e.g., Crop, Urban Industrial, Hard Linear) may be appropriate.
#' 
#' 
#' @param model Vector defined as "Vegetation", "Soil", or "All" so the appropriate model is adjusted.
#' @param coef Vector of characters that the user wishes to modify.
#' @param value User defined value to be assigned to the specified coefficients. Adjustment is made across all boostraps.
#'
#' @import Matrix
#'
#' @export
#'

coefficient_adjustment <- function(model, coef = NULL, value = NULL) {
    
    # Create a new object that contains the modified coefficients hidden to the user
    species.coefs.mod <- species.coefs
    
    # Loop through the modification for each taxonomic group
    for(taxon in names(species.coefs.mod)) {
        
        # Define the linkid
        link.fun <- link_function[[taxon]]
        
        # Vegetation adjustment
        if(model %in% c("Vegetation", "All")) {
            
            # If invalid coefficients, stop and provide warning message
            if(!all(coef %in% colnames(species.coefs.mod[[taxon]]$vegetation))) {
                
                invalid.coefs <- coef[!(coef %in% colnames(species.coefs.mod[[taxon]]$vegetation))]
                
                stop(paste("Unknown vegetation coefficient were provided:",
                           invalid.coefs,
                           collapse=", "))
                
            }

            # Adjustments
            for(coef.adjust in 1:length(coef)) {
                
                species.coefs.mod[[taxon]]$vegetation[, coef[coef.adjust], ] <- link.fun(value[coef.adjust])
                
            }
            
        }

        # Soil adjustment
        # We skip the amphibian models since none are present
        if(taxon == "Amphibians") { next }
        
        if(model %in% c("Soil", "All")) {
            
            # If invalid coefficients, stop and provide warning message
            if(!all(coef %in% colnames(species.coefs.mod[[taxon]]$soil))) {
                
                invalid.coefs <- coef[!(coef %in% colnames(species.coefs.mod[[taxon]]$soil))]
                
                stop(paste("Unknown soil coefficient were provided:",
                           invalid.coefs,
                           collapse=", "))
                
            }
            
            # Adjustments
            for(coef.adjust in 1:length(coef)) {
                
                species.coefs.mod[[taxon]]$soil[, coef[coef.adjust], ] <- link.fun(value[coef.adjust])
                
            }
            
        }
        
    }
    
    assign(x = ".species.coefs.mod", 
           value = species.coefs.mod, 
           envir = .GlobalEnv)
    
    stop("Coefficients have been modified based on user defined inputs.")
    
}

