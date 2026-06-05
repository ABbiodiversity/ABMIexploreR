#' Mammal predictions
#'
#' @description Function for calculating the species level predictions for mammals
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

.mammal_predicts <- function(species, taxon = NULL, veg = NULL, 
                                  soil = NULL, climate = NULL, boot = 0,
                                  modified = FALSE) {
    
    # If modified coefficients are defined by the user, use them
    if(modified) {
        
        COEFS <- .species.coefs.mod
        
    } else {
        
        COEFS <- species.coefs
        
    }
    
    # # If boot equals zero, pull the median bootstrap value
    # # Currently we only have a single model. Therefore, we fix boot to 1 regardless of the user defined option
    # if (boot == 0) {
    #     
    #     # Call the species lookup table function
    #     species.lookup <- abmi_species()
    #     
    #     # Isolate the relevant bootstrap iteration
    #     boot <- species.lookup[species.lookup$SpeciesID %in% species, "Bootstrap"]
    #     
    # } 
    
    boot <- 1
    
    # Define vegetation and soil coefficients
    if(species %in% rownames(COEFS[[taxon]]$Climate)) {
        
        species.climate <- COEFS[[taxon]]$Climate[species,,boot]
        
    } else {
        
        species.climate <- NULL
        
    }
    
    if(species %in% rownames(COEFS[[taxon]]$Vegetation)) {
        
        species.veg <- COEFS[[taxon]]$Vegetation[species,,boot]
        
    } else {
        
        species.veg <- NULL
        
    }
    
    if(species %in% rownames(COEFS[[taxon]]$Soil)) {
        
        species.soil <- COEFS[[taxon]]$Soil[species,,boot]
        
    } else {
        
        species.soil <- NULL
        
    }
    
    # Define link and inverse link functions
    # These are separate for the total abundance and climate
    inv_link_climate <- inv_link_function[[taxon]]$Climate
    link_climate <- link_function[[taxon]]$Climate
    
    inv_link_abundance <- inv_link_function[[taxon]]$TotalAbundance
    link_abundance <- link_function[[taxon]]$TotalAbundance
    
    # Define blank object to be returned to user
    veg.pred <- NULL
    soil.pred <- NULL
    
    #######################################
    # Create the global bioclimatic model #
    #######################################
    
    # Standardize the climate data
    climate <- climate$Climate
    climate.global <- as.matrix(climate[, names(species.climate)])
    climate.coef <- species.climate[colnames(climate.global)]
    
    # Predict space/climate component
    climate.global <- matrix(inv_link_climate(drop(climate.global %*% climate.coef)), ncol = 1,
                           dimnames = list(rownames(climate.global), "Climate"))
    
    # Truncate climate prediction
    climate.global <- ifelse(climate.global >= quantile(climate.global, 0.99),
                           quantile(climate.global, 0.99),
                           climate.global)
    
    #################################################
    # Perform the separate veg and soil predictions #
    #################################################
    
    # Identify if a vegetation prediction can be made
    if(!is.null(veg) & !is.null(species.veg)) {
        
        # Subset the climate data to the region of interest
        climate.pred <- climate.global[rownames(veg), ]
        
        # Use this to predict the joint climate contribution
        climate.pred <- (climate.pred * species.veg["Climate"])
        
        # Using these prediction, create a matrix and get the climate adjusted veg coefficients
        climate.matrix <- matrix(climate.pred, nrow = nrow(veg), ncol = ncol(veg))
        
        # Standardize the vegetation coefficients (might not be needed)
        veg.coef <- species.veg[colnames(veg)]
        
        # Prediction
        veg.coef <- t(t(climate.matrix) + veg.coef)
        veg.pred <- rowSums(veg * inv_link_abundance(veg.coef))
        
    }
    
    # Identify if soil prediction can be made
    if(!is.null(soil) & !is.null(species.soil)) {
        
        # Calculate the pAspen component
        paspen.pred <- matrix(climate[, "pAspen"], ncol = 1,
                              dimnames = list(rownames(climate), "pAspen"))
        
        # Subset the climate and pAspen data to the region of interest
        climate.pred <- climate.global[rownames(soil), ]
        paspen.pred <- paspen.pred[rownames(soil), ]
        
        # Use these to predict the joint climate contribution
        climate.pred <- (climate.pred * species.soil["Climate"])
        paspen.pred <- (paspen.pred * species.soil["pAspen"])
        climate.matrix <- matrix(climate.pred + paspen.pred, nrow = nrow(soil), ncol = ncol(soil))
        
        soil.coef <- species.soil[colnames(soil)]
        
        # Perform the joint prediction
        soil.coef <- t(t(climate.matrix) + soil.coef)
        soil.pred <- rowSums(soil * inv_link_abundance(soil.coef))
        
        # Garbage collect
        gc()
        
    }
    
    # Return the predictions to user
    return(list(Vegetation = veg.pred, 
                Soil = soil.pred))
    
}
