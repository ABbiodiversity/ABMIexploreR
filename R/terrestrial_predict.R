#' Terrestrial predictions
#'
#' @description Function for calculating the species level predictions for birds, bryophytes, lichens, soil mites, amphibians, and vascular plants
#'
#' @param species Unique Species ID. 
#' @param taxon Taxon associated with the unique Species ID.
#' @param veg Matrix of standardized vegetation information (row sum equals 1). If NULL and model is available, no prediction is generated.
#' @param soil Matrix of standardized soil information (row sum equals 1). If NULL and model is available, no prediction is generated.
#' @param climate Matrix of the bioclimatic variables that is matched to the vegetation and/or soil matrices.
#' @param boot Defines bootstrap iteration to generate the predictions from.
#'
#' @import Matrix
#'
#' @export
#'

.terrestrial_predicts <- function(species, taxon = NULL, veg = NULL, 
                                  soil = NULL, climate = NULL, boot = 0,
                                  modified = modified) {
    
    # If modified coefficients are defined by the user, use them
    if(modified) {
        
        COEFS <- .species.coefs.mod
        
    } else {
        
        COEFS <- species.coefs
        
    }
    
    # If boot equals zero, pull the median bootstrap value
    if (boot == 0) {
        
        # Call the species lookup table function
        species.lookup <- abmi_species()
        
        # Isolate the relevant bootstrap iteration
        boot <- species.lookup[species.lookup$SpeciesID %in% species, "Bootstrap"]
        
    } 
    
    # Define vegetation and soil coefficients
    if(species %in% rownames(COEFS[[taxon]]$climate)) {
        
        species.climate <- COEFS[[taxon]]$climate[species,,boot]
        
    } else {
        
        species.climate <- NULL
        
    }
    
    if(species %in% rownames(COEFS[[taxon]]$vegetation)) {
        
        species.veg <- COEFS[[taxon]]$vegetation[species,,boot]
        
    } else {
        
        species.veg <- NULL
        
    }
    
    if(species %in% rownames(COEFS[[taxon]]$soil)) {
        
        species.soil <- COEFS[[taxon]]$soil[species,,boot]
        
    } else {
        
        species.soil <- NULL
        
    }
    
    # Define link and inverse link functions
    inv_link <- inv_link_function[[taxon]]
    link <- link_function[[taxon]]
    
    # Define blank object to be returned to user
    veg.pred <- NULL
    soil.pred <- NULL
    
    #################################################
    # Perform the separate veg and soil predictions #
    #################################################
    
    # Identify if a vegetation prediction can be made
    if(!is.null(veg) & !is.null(species.veg)) {
        
        # Standardize the climate data
        climate.veg <- as.matrix(climate$vegetation[, names(species.climate)])
        climate.coef <- species.climate[colnames(climate.veg)]
        
        # Predict space/climate component
        climate.pred <- matrix(inv_link(drop(climate.veg %*% climate.coef)), ncol = 1,
                               dimnames = list(rownames(climate.veg), "Climate"))
        
        # Truncate climate prediction
        climate.pred <- ifelse(climate.pred >= quantile(climate.pred, 0.99),
                               quantile(climate.pred, 0.99),
                               climate.pred)
        
        # Use this to predict the joint climate contribution
        climate.pred <- (climate.pred * species.veg["Climate"])
        
        # Using these prediction, create a matrix and get the climate adjusted veg coefficients
        climate.matrix <- matrix(climate.pred, nrow = nrow(veg), ncol = ncol(veg))
        
        # Standardize the vegetation coefficients (might not be needed)
        veg.coef <- species.veg[colnames(veg)]
        
        # Prediction
        veg.coef <- t(t(climate.matrix) + veg.coef)
        veg.pred <- rowSums(veg * inv_link(veg.coef))
        
    }
    
    # Identify if vegetation or soil based models
    if(!is.null(soil) & !is.null(species.soil)) {
        
        # Standardize the climate data
        climate.soil <- as.matrix(climate$soil[, names(species.climate)])
        climate.coef <- species.climate[colnames(climate.soil)]
        
        # Predict space/climate component
        climate.pred <- matrix(inv_link(drop(climate.soil %*% climate.coef)), ncol = 1,
                               dimnames = list(rownames(climate.soil), "Climate"))
        paspen.pred <- matrix(climate$soil[, "pAspen"], ncol = 1,
                              dimnames = list(rownames(climate$soil), "pAspen"))
        
        # Truncate climate prediction
        climate.pred <- ifelse(climate.pred >= quantile(climate.pred, 0.99),
                               quantile(climate.pred, 0.99),
                               climate.pred)
        
        # Use these to predict the joint climate contribution
        climate.pred <- (climate.pred * species.soil["Climate"])
        paspen.pred <- (paspen.pred * species.soil["pAspen"])
        climate.matrix <- matrix(climate.pred + paspen.pred, nrow = nrow(soil), ncol = ncol(soil))
        
        soil.coef <- species.soil[colnames(soil)]
        
        # Perform the joint prediction
        soil.coef <- t(t(climate.matrix) + soil.coef)
        soil.pred <- rowSums(soil * inv_link(soil.coef))
        
        # Garbage collect
        gc()
        
    }
    
    # Return the predictions to user
    return(list(veg = veg.pred, 
                soil = soil.pred))
    
}
